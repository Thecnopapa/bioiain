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
    def __init__(self, *args, hidden_dims=[10, 10], num_classes=20, dropout=0, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout
        self._current_state = None
        self.set_mode("autoencoder")
        self.data["estimator"] = None


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


    def forward(self, x, to_latent=False, from_latent=False):
        if not from_latent:
            x = super().forward(x, submodel_name="encoder")
        if not to_latent:
            x = super().forward(x, submodel_name="decoder")
        return x


    def latent_generator(self, dataset):
        print("generating latent embeddings")
        with torch.no_grad():
            for n, item in enumerate(dataset):
                latent = self.forward(item.t, to_latent=True)
                yield latent.detach().numpy()


    def get_closest_latent(self, x, as_token=False):

        token = self._current_state.predict(x.detach().numpy().reshape(1, -1).astype(float))
        if as_token:
            return token
        latent = torch.Tensor(self._current_state.cluster_centers_[token])
        return latent



    def cluster_latent_space(self, dataset):
        from sklearn.cluster import KMeans

        algorithm = KMeans(n_clusters=20)
        with torch.no_grad():
            algorithm.fit(list(self.latent_generator(dataset)))

        self._current_state = algorithm
        self.data["estimator"] = self._current_state.get_params(deep=True)


    def plot_current_state(self, dataset=None):
        from ..visualisation.plots import fig2D
        from sklearn.decomposition import PCA
        from PIL import Image

        pca = PCA(n_components=2)
        state = pca.fit_transform(self._current_state.cluster_centers_.copy())

        print("Component:")
        for n, c in enumerate(pca.components_):
            print(f" -PC{n+1}: {c}")

        fig, ax = fig2D()

        for n, s in enumerate(state):
            ax.scatter(*s, color=f"C{n}")
            ax.text(*s, n)

        if dataset is not None:
            transformed = pca.transform(list(self.latent_generator(dataset)))

            for n, (e, l) in enumerate(zip(transformed, self._current_state.labels_)):
                #token = self.get_closest_latent(e, as_token=True)
                ax.scatter(*e, color=f"C{l}")



        os.makedirs("./latents", exist_ok=True)
        fig_path = f"./latents/{self}_E{self.data["epoch"]}.png"
        fig.savefig(fig_path)
        img = Image.open(fig_path)
        img = torchvision.transforms.v2.functional.pil_to_tensor(img)
        if self.writer is not None:
            self.writer.add_image(f"latent", img, self.data["epoch"])

    def predict_tokens(self, embeddings):
        from sklearn.cluster import KMeans

        estimator = KMeans()
        estimator.set_params(self.data["estimator"])
        preds = estimator.predict(embeddings)
        return preds














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





