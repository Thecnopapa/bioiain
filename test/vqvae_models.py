import os, json, sys, random
sys.path.append('..')

import torchvision.transforms.v2.functional


from src.bioiain.utilities.exceptions import *
from src.bioiain.utilities.sequences import *

from src.bioiain.utilities.maths import *

from src.bioiain.machine import DEVICE, tensor_to_numpy
from src.bioiain.machine.losses import *
from src.bioiain.machine.models import BaseModel
from src.bioiain.machine.layers import *

import matplotlib as mpl
import matplotlib.pyplot as plt
from src.bioiain.visualisation.plots import grid2D, fig2D

from PIL import Image
from sklearn.decomposition import PCA





class Summer(BaseModel):
    def __init__(self, *args, hidden_dims=None, num_classes=20, **kwargs):
        self.data={}
        if hidden_dims is None:
            hidden_dims = [kwargs.get("in_shape")[-1]]

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["latent_dims"] = self.data["hidden_dims"][-1]


        super().__init__(*args, **kwargs)


        self.layers["encoder"] = {
            "en_linear": nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-1])
        }
        self.layers["decoder"] = {
            "de_linear": nn.Linear(self.data["hidden_dims"][-1], self.data["in_shape"][0])
        }

        self.layers["codebook"] = {
            "codebook": Codebook(20, self.data["latent_dims"])
        }

        self.layers["autoencoder"] = {
            **self.layers["encoder"],
            **self.layers["codebook"],
            **self.layers["decoder"],
        }
        self.codebook_index = list(self.layers["autoencoder"].keys()).index("codebook")
        self.MSE = nn.MSELoss()

        self.optimisers.pop("default")


        self.optimisers["encoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["encoder"]
        }
        self.optimisers["decoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["decoder"]
        }
        self.optimisers["autoencoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["autoencoder"]
        }
        self.criterions["autoencoder"] = VQLoss()

        self.running_loss["encoder"] = 0
        self.running_loss["decoder"] = 0
        self.running_loss["autoencoder"] = 0

        self.set_mode("autoencoder")
        self.data["name"] = str(self)

    def __str__(self):
        return f"{self.__class__.__name__}{self.data['latent_dims']}D_{self.data['dataname']}"




    def _predict_from_latent(self, x):
        """
        Predict token form latent space tensor
        :param x: latent space tensor
        :return: index, score, latent(token), latent(embedding)/input
        """
        self.set_mode("codebook", quiet=True)
        z = self._forward(x)
        index = int(self.submodels["autoencoder"][self.codebook_index].last_index[0].detach().cpu().numpy())
        score = float(self.submodels["autoencoder"][self.codebook_index].last_loss.detach().cpu().numpy())
        return index, score, z, x

    def _predict(self, x) -> tuple[int, float, torch.Tensor, torch.Tensor]:
        """
        Predict token form embedding tensor
        :param x: embedding tensor
        :return: index, score, latent(token), latent(embedding)
        """
        x = self._encode(x)
        self.set_mode("codebook", quiet=True)
        z = self._forward(x)
        index = int(self.submodels["autoencoder"][self.codebook_index].last_index[0].detach().cpu().numpy())
        score = float(self.submodels["autoencoder"][self.codebook_index].last_loss.detach().cpu().numpy())
        return index, score, z, x

    def _encode(self, x):
        self.set_mode("encoder", quiet=True)
        x = x.to(DEVICE)
        z = self._forward(x)
        return z

    def _decode(self, x):
        self.set_mode("decoder", quiet=True)
        x = x.to(DEVICE)
        z = self._forward(x)
        zn = torch.clamp(z, min=0, max=1)
        return zn

    def _autoencode(self, x):

        x = x.to(DEVICE)
        y = self._encode(x)
        z = self._decode(y)

        return z, y

    def forward(self, x):
        # print("FORWARD")
        self.set_mode("autoencoder", quiet=True)

        x = x.to(DEVICE)
        z = self._forward(x)
        # zn = torch.clamp(z, min=0, max=1)

        # print(self.submodels["autoencoder"])
        encoding_loss = self.submodels["autoencoder"][self.codebook_index].last_loss
        decoding_loss = self.MSE(x, z)
        self.running_loss["encoder"] += encoding_loss.item()
        self.running_loss["decoder"] += decoding_loss.item()
        # print("encoding loss:", encoding_loss)
        # print("decoding loss:", decoding_loss)
        loss = self.loss(encoding_loss, decoding_loss)
        # print("loss:", loss)

        return loss, encoding_loss, decoding_loss

    def _latent_distance_matrix(self, normalise=True, discretise: int = None):
        log(2, "Generating latent distance matrix...")
        codebook = self.submodels["autoencoder"][self.codebook_index]
        tokens = np.array(list(zip(*codebook.codebook.weight.t().detach().cpu().numpy())))
        distances = {}
        for n1, t1 in enumerate(tokens):
            distances[n1, n1] = 0
            if n1 == len(tokens) - 1:
                break
            for n2, t2 in enumerate(tokens[n1 + 1:]):
                n2 = n1 + n2 + 1
                if n1 != n2:
                    d = multidimensional_distance(t1, t2)
                    distances[n1, n2] = d
                    # distances[n2, n1] = d

        max_dist = max(distances.values())
        if normalise:
            for k, d in distances.items():
                distances[k] = d / max_dist
            max_dist = 1

        if discretise is not None:
            bins = np.linspace(0, max_dist, discretise + 1)
            print(bins)
            print(len(bins))
            for k, d in distances.items():
                distances[k] = bins[np.digitize([d], bins[1:])][0]

        distances = {k: v for k, v in
                     sorted([(kk, vv) for kk, vv in distances.items()], key=lambda x: intto1(x[0]))}
        for k, d in distances.items():
            print(f"{k[0] + 1}-{k[1] + 1} ({intto1(k[0])}-{intto1([1])}): {distances[k]:2.1f}")
            pass
        return distances, max_dist

    def _build_blossum(self, discretisation=10):
        distances, max_dist = self._latent_distance_matrix(normalise=True, discretise=discretisation)
        blossum_folder = os.path.join(self.data["folder"], "matrixes")
        blossum_path = os.path.join(blossum_folder, f"matrix_{self}_E{self.data['epoch']}.mat")
        os.makedirs(blossum_folder, exist_ok=True)
        print("open", blossum_path)
        with open(blossum_path, "w") as f:
            letters = sorted(list(set([intto1(k[0]) for k in distances.keys()])))
            header = "   " + "  ".join(letters + ["X", "*"])
            f.write(header)
            f.write("\n")

            toks = sorted(list(set([k[0] for k in distances.keys()])), key=lambda x: intto1(x)) + ["X", "*"]

            upper = discretisation // 2
            lower = upper - discretisation

            for t1 in toks:
                l1 = intto1(t1)
                line = f"{l1}"
                for t2 in toks:
                    try:
                        if t1 == t2:
                            v = 0
                        elif t2 < t1:
                            v = (distances[(t2, t1)] * discretisation)
                        else:
                            v = (distances[(t1, t2)] * discretisation)
                    except:
                        v = discretisation

                    v = discretisation - v + lower
                    line += f" {v:2.0f}"
                f.write(line)
                f.write("\n")
        return blossum_path

    def plot_latent_dimensions(self, dataset, name="dimensioons", max_points=100, save=True, show=False,
                               fig_dir=None, plot_raw=True, r_threshold=5, only: int | list[int] = None):
        with torch.no_grad():
            log(1, "Plotting latent dimensions...")

            colorbar = mpl.colormaps["plasma"]

            embedding_size = dataset.get(1).t.size()[-1]
            latent_size = self.data["latent_dims"]
            latent_size_sqrt = math.ceil((latent_size + 1) ** 0.5)

            fig, axes = grid2D(latent_size_sqrt, latent_size_sqrt, height=latent_size_sqrt, width=latent_size_sqrt)
            last_ax = axes[latent_size]
            for ax in axes[latent_size:]:
                ax.set_axis_off()
            axes = axes[:latent_size]

            indexes = range(len(dataset))
            if len(dataset) > max_points:
                indexes = sorted(random.sample(list(indexes), max_points))

            points = [self._predict(item.t)[3] for n, item in enumerate(dataset) if n in indexes]
            decoded = [self._decode(p) for p in points]
            points = [tensor_to_numpy(p) for p in points]
            decoded = [tensor_to_numpy(d) for d in decoded]

            violin_settings = dict(
                showmeans=False,
                showmedians=False,
                showextrema=False,
                orientation="horizontal",
                positions=[0.5],
            )

            names = dataset.data["param_names"]

            for i, ax in enumerate(axes):
                ax.set_title(f"Dimension {i}", size=10)
                v = ax.violinplot([p[i] for p in points], **violin_settings)
                for b in v["bodies"]:
                    b.set_facecolor("black")
                    b.set_edgecolor("black")
                    b.set_alpha(0.1)
                # for f in range(embedding_size):
                #     vv = ax.violinplot([d[f] for d in data], **violin_settings)
                #     for bb in vv["bodies"]:
                #         bb.set_facecolor(f"C{f}")
                #         bb.set_edgecolor(f"C{f}")
                #         bb.set_alpha(0.1)

                sorted_data = sorted(list(zip(points, decoded)), key=lambda x: x[0][i])

                pp = [p[0][i] for p in sorted_data]

                for f in range(embedding_size):
                    # if i == len(axes):
                    #     legend_lines.append(mpl.lines.Line2D([0], [0], color=f"C{f}"))
                    #     legend_names.append(names[f] if f < len(names) else "")

                    if type(only) is list:
                        if not f in only:
                            continue
                    elif type(only) is int:
                        if f != only:
                            continue

                    dd = [d[1][f] for d in sorted_data]
                    a, b, c, d = (0, 0, 0, 0)
                    (d, c, b, a), r, rr, rrr, rrrr = np.polyfit(pp, dd, deg=3, full=True)
                    # print(r, rr, rrr, rrrr)

                    if r[0] <= r_threshold:
                        x_seq = np.linspace(min(pp), max(pp), 100)

                        ax.plot(x_seq, a + b * x_seq + c * (x_seq ** 2) + d * (x_seq ** 3), color=f"C{f}",
                                linewidth=3, alpha=1)
                        if plot_raw:
                            ax.plot(pp, dd, color=f"C{f}", alpha=0.2)
                        ax.set_ylim(-0.1, 1.1)
            for f in range(embedding_size):
                last_ax.plot(0, 0, alpha=1, linewidth=3, label=names[f], color=f"C{f}")
                last_ax.legend()

            plt.subplots_adjust()
            plt.tight_layout()

            if save:
                if fig_dir is None:
                    fig_dir = os.path.join(self.data["folder"], name)
                os.makedirs(fig_dir, exist_ok=True)
                fig_path = os.path.join(fig_dir,
                                        f"{name}_{self}_E{self.data['epoch']}{'_raw' if plot_raw else ''}.png")
                log(1, "Saving to: open", fig_path)
                fig.savefig(fig_path)
            if show:
                fig.show()
                plt.show(block=True)
            plt.close(fig)

            if self.writer is not None and save:
                img = Image.open(fig_path)
                img = torchvision.transforms.v2.functional.pil_to_tensor(img)
                if plot_raw:
                    self.writer.add_image(f"dimensions/raw", img, global_step=self.data["epoch"])
                else:
                    self.writer.add_image(f"dimensions/clean", img, global_step=self.data["epoch"])
                del img

    def plot_latent_space(self, dataset=None, letters=True, seed=6, fig_dir=None, show=False, plot_preds=None,
                          max_points=1000, mesh_points=None, save=True):
        with torch.no_grad():
            log(1, "Plotting latent space...")
            colorbar = mpl.colormaps["plasma"]

            mesh = False
            if mesh_points is not None:
                mesh = True

            if dataset is None and plot_preds is None:
                fig, ax = fig2D(figsize=[3000, 3000])
                axes = []
            else:
                if dataset is not None:
                    size_emb = dataset.get(1).t.size()[-1] + 1
                else:
                    size_emb = plot_preds[0][3].size()[-1] + 1
                size = math.ceil(size_emb ** 0.5)
                fig, axes = grid2D(size, size)
                ax = axes[0]
                axes = axes[1:size_emb]

            codebook = self.submodels["autoencoder"][self.codebook_index]

            o_latent = np.array(list(zip(*codebook.codebook.weight.t().detach().cpu().numpy())))
            # print("Latent:")
            # print(latent)
            # print("##")

            pca = None
            if codebook.latent_dims > 2:
                # mesh = False
                log(2, "Performing PCA on latent...")
                pca = PCA(n_components=2, random_state=seed)
                latent = pca.fit_transform(o_latent)

                log(3, "PCA components:")
                pca_text = ""
                for n, c in enumerate(pca.components_):
                    pca_text += f"PC{n + 1}: " + ", ".join([f"{x:3.2f}" for x in c]) + "\n"
                    log(4, f"PC{n + 1}: {c}")
                self.add_text("pca/components", pca_text)
                # print(latent)
            else:
                latent = o_latent

            names = ["tokens"] + dataset.data["param_names"]

            if mesh:
                x_min = min([l[0] for l in latent])
                x_max = max([l[0] for l in latent])
                y_min = min([l[1] for l in latent])
                y_max = max([l[1] for l in latent])

                x_size = (x_max - x_min) / mesh_points
                y_size = (y_max - y_min) / mesh_points

                padding = mesh_points // 10
                x_padding = padding * x_size
                y_padding = padding * y_size

                x_range = np.linspace(x_min - x_padding, x_max + x_padding, mesh_points + (padding * 2))
                y_range = np.linspace(y_min - y_padding, y_max + y_padding, mesh_points + (padding * 2))

                for x in x_range:
                    for y in y_range:
                        m = (float(x - (x_size / 2))), float((y - (y_size / 2)))

                        if pca is not None:
                            m = pca.inverse_transform(np.array(m).reshape(1, -1))[0].astype("f")
                            # print(m)
                            # print(m.dtype)

                        token, _, _, _ = self._predict_from_latent(torch.tensor(m).to(DEVICE))

                        pred = self._decode(torch.tensor(m).to(DEVICE))
                        pred = pred.detach().cpu().numpy()
                        if token is not None:
                            ax.add_patch(mpl.patches.Rectangle((x, y), x_size, y_size, color=f"C{token}"))
                        for i, axx in enumerate(axes):
                            c = colorbar(round(pred[i].item() * 255))
                            axx.add_patch(mpl.patches.Rectangle((x, y), x_size, y_size, color=c))

            if dataset is not None:

                indexes = range(len(dataset))
                if len(dataset) > max_points:
                    indexes = sorted(random.sample(list(indexes), max_points))

                log(2, f"Plotting dataset... ({len(indexes)}/{len(dataset)})")
                for n, item in enumerate(dataset):
                    log(3, f"{n + 1}/{len(dataset)}", end="\r")
                    if not n in indexes:
                        continue
                    # token = self.get_closest_latent(e, only_id=True)
                    token, _, _, point = self._predict(item.t)
                    point = point.detach().cpu().numpy()
                    if codebook.latent_dims > 2:
                        # print(point)
                        point = pca.transform(point.reshape(1, -1))[0]
                        # print(point)
                    if mesh:
                        ax.scatter(*point, color=f"C{token}", edgecolors='black')
                    else:
                        ax.scatter(*point, color=f"C{token}")
                    for i, axx in enumerate(axes):
                        c = colorbar(round(item.t[i].item() * 255))
                        if mesh:
                            axx.scatter(*point, color=c, edgecolors='black')
                        else:
                            axx.scatter(*point, color=c)
                        # if random.random() < 0.05:
                        #    axx.text(*point, f"{item.t[i].item():3.2f}", color="black")

            if plot_preds is not None:
                log(2, f"Plotting predictions... ({len(plot_preds)})")
                for n, (token, _, point, _) in enumerate(plot_preds):
                    log(3, f"{n + 1}/{len(plot_preds)}", end="\r")
                    point = point.detach().cpu().numpy()
                    if codebook.latent_dims > 2:
                        point = pca.transform(point.reshape(1, -1))[0]
                    if mesh:
                        ax.scatter(*point, color=f"C{token}", edgecolors='black', marker="s")
                    else:
                        ax.scatter(*point, color=f"C{token}", marker="s")
                    for i, axx in enumerate(axes):
                        c = colorbar(round(item.t[i].item() * 255))
                        if mesh:
                            axx.scatter(*point, color=c, edgecolors='black', marker="s")
                        else:
                            axx.scatter(*point, color=c, marker="s")

            for n, s in enumerate(latent):
                for a, axx in enumerate([ax] + axes):
                    if mesh:
                        axx.scatter(*s, color=f"C{n}", edgecolors='black')
                    else:
                        axx.scatter(*s, color=f"C{n}")
                    if letters:
                        axx.text(*s, intto1(n), backgroundcolor=f"C{n}", fontsize="xx-small")
                    else:
                        axx.text(*s, n + 1, backgroundcolor=f"C{n}", fontsize="xx-small")
                    try:
                        axx.set_title(names[a])
                    except:
                        pass

            if save:
                if fig_dir is None:
                    fig_dir = os.path.join(self.data["folder"], "latents")
                os.makedirs(fig_dir, exist_ok=True)
                fig_path = os.path.join(fig_dir, f"latent_{self}_E{self.data['epoch']}.png")
                log(1, "Saving to: open", fig_path)
                fig.savefig(fig_path)
            if show:
                fig.show()
                plt.show(block=True)
            plt.close(fig)

            if self.writer is not None and save:
                img = Image.open(fig_path)
                img = torchvision.transforms.v2.functional.pil_to_tensor(img)
                self.writer.add_image(f"latent", img, global_step=self.data["epoch"])
                del img

    def plot_tokens(self):
        with torch.no_grad():
            log(1, "Plotting Tokens...")

            fig, axes = grid2D(5, 4)

            tokens = self.submodels["autoencoder"][
                self.codebook_index].codebook.weight.detach().cpu()  # .numpy().tolist()
            # print("TOKENS")
            # print(len(tokens))
            # print(tokens)
            for n, (ax, token) in enumerate(zip(axes, tokens)):
                # print(token)
                dec = self._decode(torch.Tensor(token)).detach().cpu().numpy()
                i_length, j_length, i_j_angle, i_j_length = dec[0], dec[1], dec[2], dec[3]

                i_length = i_length * 2.4
                j_length = j_length * 2.4
                i_j_angle = i_j_angle * 180
                i_j_length = i_j_length * 10

                cv1_start = (0, 0)

                cv1_end = (0, i_length)

                cv2_start = (i_j_length, 0)

                cv2_end = rotate2D(cv2_start, (i_j_length, j_length), i_j_angle)

                ax.scatter(*cv1_start)
                ax.scatter(*cv1_end)
                ax.plot(*zip(cv1_start, cv1_end))
                ax.text(*cv1_start, "I")

                ax.text(i_j_length / 2, 0, f"{i_j_angle:3.1f}°")

                ax.scatter(*cv2_start)
                ax.scatter(*cv2_end)
                ax.plot(*zip(cv2_start, cv2_end))
                ax.text(*cv2_start, "J")

                ax.set_title(
                    f"Token {n}: i:{i_length:3.2f} j:{j_length:3.2f} d:{i_j_length:3.1f} a:{i_j_angle:3.1f}°")

            save_path = os.path.join(self.data['folder'], "tokens", f"tokens_{self}_E{self.data['epoch']}.png")
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            fig.savefig(save_path)
            plt.close(fig)

            if self.writer is not None:
                img = Image.open(save_path)
                img = torchvision.transforms.v2.functional.pil_to_tensor(img)
                self.writer.add_image(f"tokens", img, global_step=self.data["epoch"])
                del img

    def _tokenise(self, dataset):
        log(1, "Generating tokens...")
        log(2, "Dataset:", dataset)
        tok_fasta_path = dataset.data["fasta_path"].replace(".fasta", f".{self.get_fname()}.tokens.fasta")

        with open(tok_fasta_path, mode="w") as f:
            for ns, (name, strucc) in enumerate(dataset.embeddings.items()):
                print(f"{ns:4d}/{dataset.n_ids():4d} {name}", end="\r")
                start_n = strucc["start"]
                end_n = strucc["end"]
                tok_seq = ""
                for n in range(start_n, end_n):
                    item = dataset.get(n)
                    pred = self._predict(item.t)
                    token_n = pred[0]
                    token = intto1(token_n)
                    tok_seq += token

                if len(tok_seq) != len(strucc["sequence"]):
                    # log("warning", f"N ({len(tok_seq)}) of generated tokens does not mach embedding sequence length ({len(strucc['sequence'])})")
                    pass
                # log(2, tok_seq)
                f.write(f"> {name}_tokens\n")
                f.write(f"{tok_seq}\n")
        log(2, "Token fasta path:", tok_fasta_path)
        if dataset.data.get("tokenised", None) is None:
            dataset.data["tokenised"] = {}
        dataset.data["tokenised"][self.get_fname()] = tok_fasta_path
        dataset.save()
        return tok_fasta_path

    def _align_tokens(self, dataset, token_fasta_path, matrix="ID", matrix_path=None, **kwargs):
        log("header", "Aligning tokens...")
        log(1, "Token fasta path:", token_fasta_path)
        log(2, "Dataset:", dataset)
        log(2, "Matrix:", matrix, f"({matrix_path})" if matrix_path else "")

        msa = CLUSTAL(token_fasta_path, verbose=True, out_folder=dataset.data["folder"], matrix=matrix,
                      matrix_path=matrix_path, build_tree=True, **kwargs)

        return msa.msa_path






