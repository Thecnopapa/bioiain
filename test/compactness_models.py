import os, json, sys, random
sys.path.append('..')

import torchvision.transforms.v2.functional


from src.bioiain.utilities.exceptions import *
from src.bioiain.utilities.sequences import *

from src.bioiain.utilities.maths import *

from src.bioiain.machine import DEVICE, tensor_to_numpy, Embedding, humanise
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


class Structure3DEmbedding(Embedding):
    iter_dim = None
    def __init__(self, img_size, subfolder=None, **kwargs):
        if subfolder is None:
            subfolder = ""
        else:
            subfolder = f"{subfolder}_"
        super().__init__(subfolder=f"{subfolder}size_{img_size}", **kwargs)

class SaProt3DEmbedding(Structure3DEmbedding):
    def __init__(self, *args, saprot_model=None, **kwargs):
        super().__init__(*args, subfolder=saprot_model, **kwargs)

class Compactness3DEembedding(Structure3DEmbedding):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, subfolder="compactness", **kwargs)

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




class Saprot3Dto1(BaseModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        in_channels = self.data["in_shape"][0]
        img_size = self.data["img_size"] = self.data["in_shape"][-1]
        hc = [in_channels*3, in_channels*2]
        self.data["hidden_channels"] = hc
        self.data["size_reduction"] =  [(4 + 2), (2 + 2)]
        triplet_size:int = img_size - self.data["size_reduction"][0]
        max_size = self.data["max_size"] = img_size - sum(self.data["size_reduction"])


        log(2, "In shape", self.data["in_shape"])
        log(2, "In channels", in_channels)
        log(2, "Hidden channels", hc)
        log(2, "Max size", humanise(max_size))


        self.layers["convolution_common"] = {
            "conv3dMixed": nn.Conv3d(
                in_channels= in_channels,
                out_channels= hc[0],
                kernel_size=4,
                stride=1,
            ),
            "conv_reluMixed": nn.ReLU(),
            "pool3dC": nn.MaxPool3d(
                kernel_size=4,
                stride=1,
            ),
        }
        self.layers["convolution_mixed"] = {
            "conv3dM": nn.Conv3d(
                in_channels=hc[0],
                out_channels=hc[0]*2,
                kernel_size=2,
                stride=1,
            ),
            "conv_reluM": nn.ReLU(),
            "pool3dM": nn.MaxPool3d(
                kernel_size=4,
                stride=1,
            ),
        }
        self.layers["convolution_x"] = {
            "conv3dX": nn.Conv3d(
                in_channels=in_channels,
                out_channels=hc[1],
                #kernel_size=[2, 2, triplet_size],
                kernel_size=2,
                stride=1,
            ),
            "conv_reluX": nn.ReLU(),
            "pool3dX": nn.MaxPool3d(
                #kernel_size=[4, 4, 1],
                kernel_size=4,
                stride=1,
            ),
        }
        self.layers["convolution_y"] = {
            "conv3dY": nn.Conv3d(
                in_channels=in_channels,
                out_channels=hc[1],
                #kernel_size=[2, triplet_size, 2],
                kernel_size=2,
                stride=1,
            ),
            "conv_reluY": nn.ReLU(),
            "pool3dY": nn.MaxPool3d(
                #kernel_size=[4, 1, 4],
                kernel_size=4,
                stride=1,
            ),
        }
        self.layers["convolution_z"] = {
            "conv3dZ": nn.Conv3d(
                in_channels=in_channels,
                out_channels=hc[1],
                #kernel_size=[triplet_size, 2, 2],
                kernel_size=2,
                stride=1,
            ),
            "conv_reluZ": nn.ReLU(),
            "pool3dZ": nn.MaxPool3d(
                #kernel_size=[1, 4, 4],
                kernel_size=4,
                stride=1,
            ),
        }

        self.layers["linear"] = {
            "linear_relu1": nn.ReLU(),
            "flatten1": nn.Flatten(),
            "conv1d1": nn.Conv1d(hc[0] * 2, in_channels, kernel_size=1, stride=1),
            "linear_relu2": nn.ReLU(),
            "conv1d2": nn.Conv1d(in_channels, 1, kernel_size=1, stride=1),
        }

        self.layers["classifier"] = {
            "lineal_classifier": nn.Linear(max_size**3, 1),
        }

        self.layers["decoder"] = {

        }
        self.layers["compressor"] = {
            "compressor1": nn.MaxPool3d(
                kernel_size=4,
                stride=2,
            ),
            "compressor2": nn.MaxPool3d(
                 kernel_size=2,
                 stride=1,
             ),
            # "compressor3": nn.MaxPool3d(
            #     kernel_size=2,
            #     stride=1,
            # ),
            # "compressor4": nn.MaxPool3d(
            #     kernel_size=4,
            #     stride=1,
            # ),
        }


        self.layers["default"] = {
            **self.layers["convolution_common"], # --> [3840, 10, 10, 10]
            **self.layers["convolution_mixed"],
            #**self.layers["convolution_x"], # --> [2560, 6, 6, 6] x3
            #**self.layers["convolution_y"],
            #**self.layers["convolution_z"], # --> [7680, 6, 6, 6] (total)
            **self.layers["linear"], # --> [1, 216]
            **self.layers["classifier"],
        }

        self.criterions["default"] = ImgAndClassifierLoss()

    def forward(self, x, classifier=True):
        self.set_mode("default")
        x = self._forward(x, "convolution_common")
        if self.mode == "splitted":
            x, y, z = torch.split(x, self.data["in_shape"][0])
            x, y, z = self._forward(x, "convolution_x"), self._forward(y, "convolution_y"), self._forward(z, "convolution_z")
            # print(x.shape, y.shape, z.shape)
            # x, y, z = x.reshape(x.shape[0], x.shape[1]*x.shape[2]*x.shape[3]), y.reshape(y.shape[0], y.shape[1]*y.shape[2]*y.shape[3]), z.reshape(z.shape[0], z.shape[1]*z.shape[2]*z.shape[3])
            # print(x.shape, y.shape, z.shape)
            x = torch.cat((x, y, z))
        else:
            x = self._forward(x, "convolution_mixed")

        #print(x.shape)
        x = self._forward(x, "linear")
        if classifier:
            c = self._forward(x, "classifier")
            c = c.reshape(1)
        else:
            c = None
        x = x.reshape(1, self.data["max_size"], self.data["max_size"], self.data["max_size"])
        return x, c

    def compress(self, x):
        with torch.no_grad():
            #print("compressing:", x.shape)
            x = self._forward(x, "compressor")
            #print("compressed:", x.shape)
            return x
