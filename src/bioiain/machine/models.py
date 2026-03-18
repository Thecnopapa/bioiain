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

from .losses import *
from .base_model import BaseModel
from .base_model import BaseModel as CustomModel # Compatibility





class Despair(BaseModel):
    def __init__(self, *args, hidden_dims=[50, 10], num_classes=20, dropout=0, **kwargs):
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


    def latent_generator(self, dataset):
        log(2, "Generating latent embeddings...")
        n_items = len(dataset)
        for n, item in enumerate(dataset):
            print(f"{n}/{n_items}", end="\r")
            latent = self.forward(item.t, to_latent=True)
            yield latent.detach().cpu().numpy()


    def get_closest_latent(self, x, as_token=False):

        token = self._current_state.predict(x.detach().cpu().numpy().reshape(1, -1).astype(float))
        score = self._current_state.score(x.detach().cpu().numpy().reshape(1, -1).astype(float))
        if as_token:
            return token, score
        latent = torch.tensor(self._current_state.cluster_centers_[token], requires_grad=True).float().to(DEVICE)
        return latent, score



    def cluster_latent_space(self, dataset, seed=6):
        log(1, "CLustering latent space...")
        from sklearn.cluster import KMeans
        import pickle

        algorithm = KMeans(n_clusters=20, random_state=seed)
        with torch.no_grad():
            algorithm.fit(list(self.latent_generator(dataset)))

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

        pca = PCA(n_components=2, random_state=seed)
        state = pca.fit_transform(self._current_state.cluster_centers_.copy())

        log(2, "PCA components:")
        for n, c in enumerate(pca.components_):
            log(3, f"PC{n+1}: {c}")

        fig, ax = fig2D()

        for n, s in enumerate(state):
            ax.scatter(*s, color=f"C{n}")
            ax.text(*s, n)

        if dataset is not None:
            with torch.no_grad():
                transformed = pca.transform(list(self.latent_generator(dataset)))

            for n, (e, l) in enumerate(zip(transformed, self._current_state.labels_)):
                #token = self.get_closest_latent(e, as_token=True)
                ax.scatter(*e, color=f"C{l}")


        fig_dir = os.path.join(self.data["folder"], "latents")
        os.makedirs(fig_dir, exist_ok=True)
        fig_path = os.path.join(fig_dir, f"{self}_E{self.data["epoch"]}.png")
        fig.savefig(fig_path)
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
        from ..visualisation.plots import grid2D
        from PIL import Image

        fig, axes = grid2D(5, 4)

        print(axes)


        for n, ax in enumerate(axes):
            print(n, ax)
            self.draw_token(n, ax=ax)

        if show:
            fig.show()
        if save:
            save_path = os.path.join(self.data["folder"], "tokens", f"tokens.png")
            os.makedirs(os.path.dirname(save_path), exist_ok=True)
            fig.savefig(save_path)

            if self.writer is not None:
                img = Image.open(save_path)
                img = torchvision.transforms.v2.functional.pil_to_tensor(img)
                self.writer.add_image(f"tokens", img, global_step=self.data["epoch"])
                del img

    def draw_token(self, token_id, show=False, save=False, fig=None, ax=None):
        from ..visualisation.plots import fig2D
        from ..utilities.maths import rotate2D
        print("Drawing token:", token_id)
        print("ax:", ax)


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
        ax.text(*cv1_start, "i")

        ax.text(i_j_length/2, 0, f"{i_j_length:3.1}")

        ax.scatter(*cv2_start)
        ax.scatter(*cv2_end)
        ax.plot(*zip(cv2_start, cv2_end))
        ax.text(*cv1_start, "j")

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





