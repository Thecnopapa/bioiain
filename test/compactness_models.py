import os, json, sys, random
sys.path.append('..')

import torchvision.transforms.v2.functional


from src.bioiain.utilities.exceptions import *
from src.bioiain.utilities.sequences import *

from src.bioiain.utilities.maths import *

from src.bioiain.machine import DEVICE, tensor_to_numpy, Embedding
from src.bioiain.machine.losses import *
from src.bioiain.machine.models import BaseModel
from src.bioiain.machine.layers import *

import matplotlib as mpl
import matplotlib.pyplot as plt
from src.bioiain.visualisation.plots import grid2D, fig2D

from PIL import Image
from sklearn.decomposition import PCA





class SaProtEmbedding(Embedding):
    iter_dim = 1
    pass


class structure3DEmbedding(Embedding):
    iter_dim = 0
    def __init__(self, img_size, **kwargs):
        super().__init__(subfolder=f"size_{img_size}", **kwargs)

class SaProt3DEmbedding(structure3DEmbedding):
    pass

class compactness3Dembedding(SaProt3DEmbedding):
    pass

class CompactnessMLPmk1(BaseModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        n = self.data["in_shape"][0]
        self.data["hidden_dims"] = [n*2, n*2, n]
        self.layers["default"] = {
            "linear1": nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-3]),
            "en_relu1": nn.ReLU(),
            "linear2": nn.Linear(self.data["hidden_dims"][-3], self.data["hidden_dims"][-2]),
            "en_relu2": nn.ReLU(),
            "linear3": nn.Linear(self.data["hidden_dims"][-2], self.data["hidden_dims"][-1]),
            "en_relu3": nn.ReLU(),
            "linear4": nn.Linear(self.data["hidden_dims"][-1], 1),
        }

class Compactness3Dmk1(BaseModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        n = self.data["in_shape"][0]
        self.data["hidden_dims"] = [n, n*4, n*8, n]
        cn = self.data["hidden_dims"]
        max_size = (self.data["in_shape"][1] - 4)**3
        print(max_size)

        self.layers["convolution"] = {
            "conv3d1": nn.Conv3d(
                in_channels= cn[0],
                out_channels= cn[1],
                kernel_size=4,
                stride=1,
            ),
            "conv_relu1": nn.ReLU(),
            "conv3d2": nn.Conv3d(
                in_channels=cn[1],
                out_channels=cn[2],
                kernel_size=2,
                stride=1,
            )
        }
        self.layers["linear"] = {
            "flatten1": nn.Flatten(),
            "linear_relu1": nn.ReLU(),
            "conv1d": nn.Conv1d(8, 1, kernel_size=1, stride=1),
            "linear_relu2": nn.ReLU(),
            "linear": nn.Linear(max_size, 1280),
            "linear_relu3": nn.ReLU(),
        }
        self.layers["classifier"] = {
            "classifier_head": nn.Linear(1280, 1),
        }
        self.layers["default"] = {
            **self.layers["convolution"],
            **self.layers["linear"],
            **self.layers["classifier"],
        }
