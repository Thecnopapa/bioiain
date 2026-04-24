import os, json, platform

import torchvision.transforms.v2.functional
from sklearn.metrics import confusion_matrix
import PIL

from ..utilities.logging import *
from ..utilities.parallel import *

import torch
import torch.nn as nn

import pandas as pd
import numpy as np

from .datasets import Item, Dataset

from torch.utils.tensorboard import SummaryWriter
import datetime
import psutil
from . import DEVICE

from ..utilities.exceptions import *
from ..utilities.sequences import *

from .losses import *
from .base_model import BaseModel
from .base_model import BaseModel as CustomModel # Compatibility
from .layers import *

import matplotlib.pyplot as plt




class Despair(BaseModel):
    def __init__(self, *args, hidden_dims=[10, 2], num_classes=20, dropout=0, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout
        self._current_state = None
        self.set_mode("autoencoder")
        self.data["tokens"] = None


        self.layers["encoder"] = {
            "en_l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "en_relu1": nn.ReLU(),
            "en_l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
        }

        self.layers["decoder"] = {
            "de_l2": nn.Linear(hidden_dims[1], hidden_dims[0]),
            "de_relu1": nn.ReLU(),
            "de_l1": nn.Linear(hidden_dims[0], self.data["in_shape"][0]),


        }

        self.optimisers.pop("default")
        self.optimisers["autoencoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["encoder", "decoder"]
        }
        self.optimisers["encoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["encoder"]
        }
        self.optimisers["decoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["decoder"]
        }

        self.criterions["autoencoder"] = VQLoss()
        self.running_loss["encoder"] = 0
        self.running_loss["decoder"] = 0


    def forward(self, x, to_latent=False, from_latent=False):
        x = x.to(DEVICE)
        if not from_latent:
            x = super().forward(x, submodel_name="encoder")
        if not to_latent:
            x = super().forward(x, submodel_name="decoder")
        return x


    def train(self, item, i=None, n_items=None):

        item.to(DEVICE)
        latent = self(item.t, to_latent=True)
        #print("latent", latent)
        latent, token_latent, origin_latent = self.get_closest_latent(latent)
        #print("token latent", token_latent)

        encoder_loss = self.criterions["autoencoder"].encoder_loss(latent, token_latent, origin_latent, commitment=0)
        self.running_loss["encoder"] += encoder_loss.item()

        #print("encoder_loss", encoder_loss)

        latent_diff = torch.sub(latent, token_latent)
        #print("latent diff", latent_diff)
        new_latent = torch.sub(latent, latent_diff)
        #print("new latent", new_latent)
        out = self(new_latent.to(DEVICE), from_latent=True)

        decoder_loss = self.criterions["autoencoder"].decoder_loss(out, item.t)
        self.running_loss["decoder"] += decoder_loss.item()


        loss = self.loss(encoder_loss, decoder_loss)

        if i is not None and n_items is not None:

            print(f"{i}/{n_items} LOSS: {loss.item():7.3f} ({encoder_loss.item():7.3f}/{decoder_loss.item():7.3f}) av:{self.running_loss[self.mode]/self.running_loss["total"]:7.3f}", end="\r")


        return loss





    def latent_generator(self, dataset):
        log(2, "Generating latent embeddings...")
        n_items = len(dataset)
        for n, item in enumerate(dataset):
            print(f"{n+1}/{n_items}", end="\r")
            latent = self.forward(item.t, to_latent=True)
            yield latent.detach().cpu().numpy()


    def get_closest_latent(self, x, only_id=False):
        from ..utilities.maths import multidimensional_com, multidimensional_distance

        target = x.detach().cpu().numpy().reshape(1, -1).astype(float)
        token_id = self._current_state.predict(target)

        token_latent = self._current_state.cluster_centers_[token_id][0]
        #print("TOKEN LATENT" ,token_latent)

        origin = multidimensional_com(self._current_state.cluster_centers_)

        distance_from_origin = multidimensional_distance(origin, target)
        distance_from_token = multidimensional_distance(target, token_latent)



        if only_id:
            return token_id, distance_from_token, distance_from_origin

        latent = torch.tensor(token_latent, requires_grad=True).float().to(DEVICE)
        origin = torch.tensor(origin, requires_grad=True).float().to(DEVICE)

        return x, latent , origin



    def cluster_latent_space(self, dataset, seed=6):
        log(1, "CLustering latent space...")
        from sklearn.cluster import KMeans
        import pickle

        algorithm = KMeans(n_clusters=20, random_state=seed)
        with torch.no_grad():
            a = np.array(list(self.latent_generator(dataset))).astype(float)
            log(2, "Fitting...")
            algorithm.fit(a)
            tracemalloc_top()
            del a

        self._current_state = algorithm
        self.data["tokens"] = self._current_state.cluster_centers_.copy().tolist()

        estimator_path = os.path.join(self.data["folder"], self.get_fname()+ ".estimator.pkl")
        pickle.dump(algorithm, open(estimator_path, "wb"))
        self.data["estimator_path"] = estimator_path

    def plot_current_state(self, dataset=None, seed=6, tokens=True):
        log(1, "Plotting current state...")
        from ..visualisation.plots import fig2D
        from sklearn.decomposition import PCA
        from PIL import Image



        fig, ax = fig2D()


        if self.data["hidden_dims"][-1] > 2:
            pca = PCA(n_components=2, random_state=seed)
            state = pca.fit_transform(self._current_state.cluster_centers_.copy())

            log(2, "PCA components:")
            for n, c in enumerate(pca.components_):
                log(3, f"PC{n+1}: {c}")

        else:
            state = self._current_state.cluster_centers_.copy()

        for n, s in enumerate(state):
            ax.scatter(*s, color=f"C{n}")
            ax.text(*s, n)

        if dataset is not None:
            with torch.no_grad():
                if self.data["hidden_dims"][-1] > 2:
                    transformed = pca.transform(list(self.latent_generator(dataset)))
                else:
                    transformed = self.latent_generator(dataset)

            for n, (e, l) in enumerate(zip(transformed, self._current_state.labels_)):
                #token = self.get_closest_latent(e, only_id=True)
                ax.scatter(*e, color=f"C{l}")


        fig_dir = os.path.join(self.data["folder"], "latents")
        os.makedirs(fig_dir, exist_ok=True)
        fig_path = os.path.join(fig_dir, f"latent_{self}_E{self.data["epoch"]}.png")
        fig.savefig(fig_path)
        plt.close(fig)
        print("saving to:", fig_path)

        if self.writer is not None:
            img = Image.open(fig_path)
            img = torchvision.transforms.v2.functional.pil_to_tensor(img)
            self.writer.add_image(f"latent", img, global_step=self.data["epoch"])
            del img
        if tokens:
            self.draw_all_tokens()

    def predict_tokens(self, latents):
        preds = self.estimator.predict(latents)
        return preds

    def estimator(self):
        import pickle
        estimator = pickle.load(open(self.data["estimator_path"], "rb"))
        return estimator



    def draw_all_tokens(self, show=False, save=True):
        log(1,"Drawing all tokens...")
        from ..visualisation.plots import grid2D
        from PIL import Image

        fig, axes = grid2D(5, 4)

        #print(axes)


        for n, ax in enumerate(axes):
            #print(n, ax)
            self.draw_token(n, ax=ax)

        if show:
            fig.show()
        if save:
            save_path = os.path.join(self.data["folder"], "tokens", f"tokens_{self}_E{self.data["epoch"]}.png")
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            fig.savefig(save_path)
            plt.close(fig)

            if self.writer is not None:
                img = Image.open(save_path)
                img = torchvision.transforms.v2.functional.pil_to_tensor(img)
                self.writer.add_image(f"tokens", img, global_step=self.data["epoch"])
                del img

    def draw_token(self, token_id, show=False, save=False, fig=None, ax=None):
        from ..visualisation.plots import fig2D
        from ..utilities.maths import rotate2D
        log(2,"Drawing token:", token_id, end="\r")
        #print("ax:", ax)


        tokens = self.estimator().cluster_centers_
        token = torch.Tensor(tokens[token_id])

        with torch.no_grad():
            reconstructed = self(token, from_latent=True)
            reconstructed = [i.item() for i in reconstructed]
            i_length, j_length, i_j_angle, i_j_length = reconstructed

        i_length = i_length *2.4
        j_length = j_length *2.4
        i_j_angle = i_j_angle*180
        i_j_length = i_j_length*10


        cv1_start = (0,0)

        cv1_end = (0, i_length)

        cv2_start = (i_j_length, 0)

        cv2_end = rotate2D(cv2_start, (i_j_length, j_length), i_j_angle)




        if ax is None:
            fig, ax = fig2D()
        ax.scatter(*cv1_start)
        ax.scatter(*cv1_end)
        ax.plot(*zip(cv1_start, cv1_end))
        ax.text(*cv1_start, "I")

        ax.text(i_j_length/2, 0, f"{i_j_length:3.1}")

        ax.scatter(*cv2_start)
        ax.scatter(*cv2_end)
        ax.plot(*zip(cv2_start, cv2_end))
        ax.text(*cv2_start, "J")

        #ax.quiver(cv1_start, cv1_end)
        #ax.quiver(cv2_start, cv2_end)

        ax.set_title(f"Token: i:{i_length:3.2f} j:{j_length:3.2f} d:{i_j_length:3.1f} a:{i_j_angle:3.1f}")


        if fig is not None:
            if show:
                fig.show()
            if save:
                save_path = os.path.join(self.data["folder"], "tokens", f"token_{token_id}.png")
                os.makedirs(os.path.dirname(save_path), exist_ok=True)

                fig.savefig(save_path)
            plt.close(fig)



class DespairLess(Despair):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.layers["encoder"] = {
            "en_linear": nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-1])
        }
        self.layers["decoder"] = {
            "de_linear": nn.Linear(self.data["hidden_dims"][-1], self.data["in_shape"][0])
        }



class Hope(DespairLess):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["codebook"] = {
            "codebook": Codebook(20, self.data["hidden_dims"][-1])
        }

        self.layers["autoencoder"] = {
            **self.layers["encoder"],
            **self.layers["codebook"],
            **self.layers["decoder"],
        }
        self.codebook_index = list(self.layers["autoencoder"].keys()).index("codebook")
        self.MSE = nn.MSELoss()

        self.optimisers["autoencoder"] = {
            "class": torch.optim.Adam,
            "layer_set": ["autoencoder"]
        }

    def _predict_from_latent(self, x):
        self.set_mode("codebook", quiet=True)
        z = self._forward(x)
        index = int(self.submodels["autoencoder"][self.codebook_index].last_index[0].detach().cpu().numpy())
        score = float(self.submodels["autoencoder"][self.codebook_index].last_loss.detach().cpu().numpy())
        return index, score, z, x


    def _predict(self, x):
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

    def forward(self, x):
        #print("FORWARD")
        self.set_mode("autoencoder", quiet=True)

        x = x.to(DEVICE)
        z = self._forward(x)
        #zn = torch.clamp(z, min=0, max=1)

        #print(self.submodels["autoencoder"])
        encoding_loss = self.submodels["autoencoder"][self.codebook_index].last_loss
        decoding_loss = self.MSE(x, z)
        self.running_loss["encoder"] += encoding_loss.item()
        self.running_loss["decoder"] += decoding_loss.item()
        #print("encoding loss:", encoding_loss)
        #print("decoding loss:", decoding_loss)
        loss = self.loss(encoding_loss, decoding_loss)
        #print("loss:", loss)

        return loss, encoding_loss, decoding_loss

    def _latent_distance_matrix(self):
        log(2, "Generating latent distance matrix...")
        codebook = self.submodels["autoencoder"][self.codebook_index]
        tokens = np.array(list(zip(*codebook.codebook.weight.t().detach().cpu().numpy())))
        distances = {}
        for n1, t1 in enumerate(tokens):
            if n1 == len(tokens) -1:
                break
            for n2, t2 in enumerate(tokens[n1+1:]):
                n2 = n1 + n2 + 1
                if n1 == n2:
                    distances[n1, n2] = 0
                else:
                    d = multidimensional_distance(t1, t2)
                    distances[n1, n2] = d
                    distances[n2, n1] = d
        print(json.dumps(distances, indent=4))
        return distances



    def plot_latent_space(self, dataset=None, seed=6, fig_dir=None, show=False, plot_preds=None, max_points=1000, mesh_points=None):
        with torch.no_grad():
            log(1, "Plotting current state...")
            from ..visualisation.plots import fig2D, grid2D
            from sklearn.decomposition import PCA
            from PIL import Image
            import matplotlib as mpl
            colorbar = mpl.colormaps["plasma"]

            mesh = False
            if mesh_points is not None:
                mesh = True

            if dataset is None:
                fig, ax = fig2D()
            else:
                size_emb = dataset.get(1).t.size()[-1]+1
                size = math.ceil(size_emb ** 0.5)
                fig, axes = grid2D(size, size)
                ax = axes[0]
                axes = axes[1:size_emb]

            codebook = self.submodels["autoencoder"][self.codebook_index]

            o_latent = np.array(list(zip(*codebook.codebook.weight.t().detach().cpu().numpy())))
            #print("Latent:")
            #print(latent)
            #print("##")

            pca = None
            if codebook.latent_dims > 2:
                mesh = False
                log(2, "Performing PCA on latent...")
                pca = PCA(n_components=2, random_state=seed)
                latent = pca.fit_transform(o_latent)

                log(3, "PCA components:")
                for n, c in enumerate(pca.components_):
                    log(4, f"PC{n+1}: {c}")

                #print(latent)
            else:
                latent = o_latent

            names = ["tokens", "len i", "len j", "angle ij", "dist ij", "dist lig", "contactability"]

            if mesh:
                x_min = min([l[0] for l in latent])
                x_max = max([l[0] for l in latent])
                y_min = min([l[1] for l in latent])
                y_max = max([l[1] for l in latent])

                x_size = (x_max - x_min) / mesh_points
                y_size = (y_max - y_min) / mesh_points

                padding = mesh_points // 10
                x_padding = padding*x_size
                y_padding = padding*y_size

                x_range = np.linspace(x_min-x_padding, x_max+x_padding, mesh_points+(padding*2))
                y_range = np.linspace(y_min-y_padding, y_max+y_padding, mesh_points+(padding*2))

                for x in x_range:
                    for y in y_range:
                        m = (float(x-(x_size/2))), float((y-(y_size/2)))
                        token = None
                        if pca is None:
                            token, _, _, _ = self._predict_from_latent(torch.tensor(m).to(DEVICE))

                        pred = self._decode(torch.tensor(m).to(DEVICE))
                        pred = pred.detach().cpu().numpy()
                        if token is not None:
                            ax.add_patch(mpl.patches.Rectangle((x,y), x_size, y_size, color=f"C{token}"))
                        for i, axx in enumerate(axes):
                            c = colorbar(round(pred[i].item() * 255))
                            axx.add_patch(mpl.patches.Rectangle((x,y), x_size, y_size, color=c))






            if dataset is not None:

                import random

                indexes = range(len(dataset))
                if len(dataset) > max_points:
                    indexes = sorted(random.sample(list(indexes), max_points))


                log(2, f"Plotting dataset... ({len(indexes)}/{len(dataset)})")
                for n, item in enumerate(dataset):
                    log(3, f"{n + 1}/{len(dataset)}", end="\r")
                    if not n in indexes:
                        continue
                    #token = self.get_closest_latent(e, only_id=True)
                    token, _, _, point = self._predict(item.t)
                    point = point.detach().cpu().numpy()
                    if codebook.latent_dims > 2:
                        #print(point)
                        point = pca.transform(point.reshape(1, -1))[0]
                        #print(point)
                    if mesh:
                        ax.scatter(*point, color=f"C{token}", edgecolors='black')
                    else:
                        ax.scatter(*point, color=f"C{token}")
                    for i, axx in enumerate(axes):
                        c = colorbar(round(item.t[i].item()*255))
                        if mesh:
                            axx.scatter(*point, color=c, edgecolors='black')
                        else:
                            axx.scatter(*point, color=c)
                        #if random.random() < 0.05:
                        #    axx.text(*point, f"{item.t[i].item():3.2f}", color="black")


            if plot_preds is not None:
                log(2, "Plotting predictions..." )
                for token, _, _, point in plot_preds:
                    ax.scatter(*point.detach().cpu().numpy(), color=f"C{token}")

            for n, s in enumerate(latent):
                for a, axx in enumerate([ax]+axes):
                    if mesh:
                        axx.scatter(*s, color=f"C{n}", edgecolors='black')
                    else:
                        axx.scatter(*s, color=f"C{n}")
                    axx.text(*s, n, backgroundcolor=f"C{n}", fontsize="xx-small")
                    try:
                        axx.set_title(names[a])
                    except:
                        pass


            if fig_dir is None:
                fig_dir = os.path.join(self.data["folder"], "latents")
            os.makedirs(fig_dir, exist_ok=True)
            fig_path = os.path.join(fig_dir, f"latent_{self}_E{self.data["epoch"]}.png")
            fig.savefig(fig_path)
            if show:
                fig.show()
            plt.close(fig)
            log(1, "Saving to: open", fig_path)


            if self.writer is not None:
                img = Image.open(fig_path)
                img = torchvision.transforms.v2.functional.pil_to_tensor(img)
                self.writer.add_image(f"latent", img, global_step=self.data["epoch"])
                del img


    def plot_tokens(self):
        with torch.no_grad():
            log(1, "Plotting Tokens...")
            from ..visualisation.plots import grid2D
            from PIL import Image


            fig, axes = grid2D(5, 4)

            tokens = self.submodels["autoencoder"][self.codebook_index].codebook.weight.detach().cpu()#.numpy().tolist()
            #print("TOKENS")
            #print(len(tokens))
            #print(tokens)
            for n, (ax, token) in enumerate(zip(axes, tokens)):
                #print(token)
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

                ax.set_title(f"Token {n}: i:{i_length:3.2f} j:{j_length:3.2f} d:{i_j_length:3.1f} a:{i_j_angle:3.1f}°")

            save_path = os.path.join(self.data["folder"], "tokens", f"tokens_{self}_E{self.data["epoch"]}.png")
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
                    token_n = pred[0] - 1
                    token = intto1(token_n)
                    tok_seq += token

                if len(tok_seq) != len(strucc["sequence"]):
                    #log("warning", f"N ({len(tok_seq)}) of generated tokens does not mach embedding sequence length ({len(strucc['sequence'])})")
                    pass
                #log(2, tok_seq)
                f.write(f"> {name}_tokens\n")
                f.write(f"{tok_seq}\n")
        log(2, "Token fasta path:", tok_fasta_path)
        if dataset.data.get("tokenised", None) is None:
            dataset.data["tokenised"] = {}
        dataset.data["tokenised"][self.get_fname()] = tok_fasta_path
        dataset.save()
        return tok_fasta_path

    def _align_tokens(self, dataset, token_fasta_path, **kwargs):
        log("header", "Aligning tokens...")
        log(1, "Token fasta path:", token_fasta_path)
        log(2, "Dataset:", dataset)

        msa = CLUSTAL(token_fasta_path, verbose=True, out_folder=dataset.data["folder"], matrix="ID", build_tree=True, **kwargs)

        return msa.msa_path








class HopeFull(Hope):
    def __init__(self, *args, hidden_dims=None, **kwargs):
        super().__init__(*args, hidden_dims=(50,50), **kwargs)


########################################################################################################################


class Golden(BaseModel):
    """
    3 Linear Model
    1280 -> 2560 -> 128
    3M params
    Custom Weighted MSE Loss
    """
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout
        self.data["weights"] = weights

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }

        self.criterions["default"] = CustomWeighted(self.data["weights"])




class GoldenTo1(Golden):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.criterions["default"] = CustomWeighted(self.data["weights"], addto1=True)


class GoldenAdamW(Golden):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.optimisers["default"]["class"] = torch.optim.AdamW




class GoldenDynamix(Golden):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.optimisers["default"]["LRS"] = customLRS

class GoldenDynamixExtra(GoldenDynamix):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.optimisers["default"]["LRS_kwargs"] = {"use_original": False}







### DEPRECATED #########################################################################################################






class DUAL_MLP_MK9(CustomModel):
    """
    MK4 with Custom LR
    """
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }


        self.optimisers["default"]["class"] = torch.optim.AdamW
        self.optimisers["default"]["LRS"] = self.customLRS
        self.optimisers["default"]["kwargs"]["fused"] = True


        self.criterions["default"] = CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])

        self._mount_submodels()

    class customLRS(torch.optim.lr_scheduler.LRScheduler):
        def __init__(self, *args, optimiser, **kwargs):
            self.o_lrs = [p["lr"] for p in optimiser.param_groups]
            super().__init__(optimiser, *args, **kwargs)

        def get_lr(self):
            print("LRS: getting lrs")
            print(self.lrs)
            return torch.Tensor(np.array(self.lrs))


        def step(self, running_loss=0.5):
            print("LRS: stepping...")
            print(running_loss)
            self.lrs = []
            for p, olr in zip(self.optimizer.param_groups, self.o_lrs):
                old_log = math.log(olr, 10)
                print("OLD_LOG", old_log)
                new_log = old_log - ((1-running_loss)*4) +2
                print("NEW_LOG", new_log)
                new_lr = 10 ** new_log
                print("NEW_LR", new_lr)
                self.lrs.append(new_lr)
                p["lr"] = torch.Tensor(np.array([new_lr]))
            print(self.lrs)
            return self.lrs




class DUAL_MLP_MK8(CustomModel):
    """
    MK4 with AdamW
    """
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }


        self.optimisers["default"]["class"] = torch.optim.AdamW
        self.optimisers["default"]["kwargs"]["fused"] = True


        self.criterions["default"] = CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])

        self._mount_submodels()





class DUAL_MLP_MK7(CustomModel):
    """
    MK6 with splitted Linear 1
    """
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.data["n_sublayers_l1"] = 10
        self.data["l_subsize_l1"] = hidden_dims[0]//self.data["n_sublayers_l1"]

        self.layers["default"] = {
            "sl1": self.splitLinear(input_dim=self.data["in_shape"][0], layer_size=self.data["l_subsize_l1"], n_layers=self.data["n_sublayers_l1"]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }

        self.criterions["default"] = CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])



        self.layers["no-dropout"] = self.layers["default"].copy()
        self.layers["no-dropout"].pop("drop1")
        self.layers["no-dropout"].pop("drop2")

        self._mount_submodels()


    class splitLinear(nn.Module):
        def __init__(self, *args, input_dim=1280, layer_size=128, n_layers=20, **kwargs):
            super().__init__(*args, **kwargs)
            self.layers = nn.ModuleList()
            self.input_dim = input_dim
            self.layer_size = layer_size
            for i in range(n_layers):
                self.layers.append(nn.Linear(self.input_dim, layer_size))


        def forward(self, x):
            #print("IN:", x.shape)
            segments = []
            for n, l in enumerate(self.layers):
                segments.append(l(x))
                #print(f"SEGMENT {n}:", segments[-1].shape)
            r = torch.cat(segments)
            #print("OUT:", r.shape)
            return r




class DUAL_MLP_MK6(CustomModel):
    """
    MK4 with no-dropout mode
    """
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }




        self.criterions["default"] = CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])



        self.layers["no-dropout"] = self.layers["default"].copy()
        self.layers["no-dropout"].pop("drop1")
        self.layers["no-dropout"].pop("drop2")

        self._mount_submodels()




class DUAL_MLP_MK5(CustomModel):
    """
    6 Lienar model (massive)
    """
    def __init__(self, *args, hidden_dims=[2560, 1280, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], hidden_dims[0]),
            "relu3": nn.LeakyReLU(),
            "l4": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu4": nn.LeakyReLU(),
            "l5": nn.Linear(hidden_dims[1], hidden_dims[0]),
            "relu5": nn.LeakyReLU(),
            "l6": nn.Linear(hidden_dims[0], hidden_dims[-1]),
            "drop2": nn.Dropout(dropout),
            "last": nn.Linear(hidden_dims[-1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }



        self.criterions["default"] = CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])

        self._mount_submodels()




class DUAL_MLP_MK4(CustomModel):
    """
    3 Linear Model
    1280 -> 2560 -> 128
    3M params
    Custom MSE Loss
    """
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }




        self.criterions["default"] = CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])

        self._mount_submodels()



class DUAL_MLP_MK3(CustomModel):
    """
    4 Linear model
    Cross Entropy Loss
    """
    def __init__(self, *args, hidden_dims=[640, 1280, 128], num_classes=4, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.ReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.ReLU(),
            "l3": nn.Linear(hidden_dims[1], hidden_dims[2]),
            "drop3": nn.Dropout(dropout),
            "relu3": nn.ReLU(),
            "l4": nn.Linear(hidden_dims[2], num_classes),
            "softmax": nn.Softmax(dim=0)
        }

        #self.criterions["default"] = self.DualLoss(self)
        self.criterions["default"] = nn.CrossEntropyLoss()

        self._mount_submodels()



class DUAL_MLP_MK2(CustomModel):
    """
    4 Linear model
    Dual Loss
    """
    def __init__(self, *args, hidden_dims=[2560, 1280, 128], num_classes=2, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)



        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.ReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.ReLU(),
            "l3": nn.Linear(hidden_dims[1], hidden_dims[2]),
            "drop3": nn.Dropout(dropout),
            "relu3": nn.ReLU(),
            "l4": nn.Linear(hidden_dims[2], num_classes),
            "hardtahn": nn.Hardtanh(min_val=0.)
        }

        self.criterions["default"] = self.DualLoss(self)

        self._mount_submodels()


    class DualLoss(object):
        def __init__(self, model):
            self.writer = model.writer
            self.model = model

        def __name__(self):
            return "DualLoss"

        def __call__(self, o, t):
            true_contact, true_outer = t[0], t[1]
            out_contact, out_outer = o[0], o[1]

            outer_loss = abs(true_outer - out_outer)
            contact_loss = abs(true_contact - out_contact)
            if "outer" not in self.model.running_loss: self.model.running_loss["outer"] = 0
            if "contactability" not in self.model.running_loss: self.model.running_loss["contactability"] = 0
            self.model.running_loss["outer"] += outer_loss
            self.model.running_loss["contactability"] += contact_loss

            return contact_loss * outer_loss



class DUAL_MLP_MK1(CustomModel):
    """
    3 Linear
    Dual Loss
    """
    def __init__(self, *args, hidden_dims=[128, 256], num_classes=2, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            #"drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            #"drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            # "softmax": nn.Softmax(dim=0)
        }

        self.criterions["default"] = self.DualLoss(self)

        self._mount_submodels()


    class DualLoss(object):
        def __init__(self, model):
            self.writer = model.writer
            self.model = model

        def __name__(self):
            return "DualLoss"

        def __call__(self, o, t):
            true_contact, true_outer = t[0], t[1]
            out_contact, out_outer = o[0], o[1]

            outer_loss = abs(true_outer - out_outer)
            contact_loss = abs(true_contact - out_contact)
            if "outer" not in self.model.running_loss: self.model.running_loss["outer"] = 0
            if "contactability" not in self.model.running_loss: self.model.running_loss["contactability"] = 0
            self.model.running_loss["outer"] += outer_loss
            self.model.running_loss["contactability"] += contact_loss

            return contact_loss + outer_loss



class MLP_MK3(CustomModel):
    def __init__(self, *args, hidden_dims=[128 ,256], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            "drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            "drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            #"softmax": nn.Softmax(dim=0)
        }


        self.criterions["default"] = self.simpleloss

        self._mount_submodels()


    @staticmethod
    def simpleloss(o, t):
            return abs(o-t)



class MLP_MK2(CustomModel):
    def __init__(self, *args, hidden_dims=[256 ,128], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            #"drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            #"drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            #"softmax": nn.Softmax(dim=0)
        }


        self.criterions["default"] = self.simpleloss

        self._mount_submodels()


    @staticmethod
    def simpleloss(o, t):
            return abs(o-t)



class MLP_MK1(CustomModel):
    def __init__(self, *args, input_dim, hidden_dims=[256 ,128], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(input_dim, hidden_dims[0]),
            "relu1": nn.ReLU(),
            "drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.ReLU(),
            "drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes)
        }

        self._mount_submodels()





