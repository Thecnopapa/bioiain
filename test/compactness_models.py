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
    def __init__(self, *args, relative=True, **kwargs):
        if relative:
            prefix=f"compactness_rel"
        else:
            prefix=f"compactness_abs"
        super().__init__(*args, subfolder=prefix, **kwargs)

class CompactnessMLPmk1(BaseModel):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        n = self.data["in_shape"][0]
        self.data["hidden_dims"] = [n*2, n*2, n]
        self.layers["default"] = {
            "linear1":    nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-3]),
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
        hc = self.data["hidden_channels"] = [in_channels*3, in_channels*2]
        latent_size = self.data["latent_size"] = 4
        conv_kernels = [4, 2]

        self.data["size_reduction"] =  [k-1 for k in conv_kernels]
        #triplet_size:int = img_size - self.data["size_reduction"][0]
        max_size = self.data["max_size"] = img_size - sum(self.data["size_reduction"])
        assert max_size >= latent_size, f"Input image is too small for this model by {latent_size - max_size}"

        first_kernel = None
        last_kernel = max(1, (max_size - latent_size) +1)
        pool_kernels = self.data["pool_kernels"] = [first_kernel, last_kernel]


        log(2, "In shape", self.data["in_shape"])
        log(2, "In channels", in_channels)
        log(2, "Hidden channels", hc)
        log(2, "Max size", humanise(max_size))
        log(2, "Pool kernels", pool_kernels)


        self.layers["convolution_common"] = {
            "conv3dMixed": nn.Conv3d(
                in_channels= in_channels,
                out_channels= hc[0],
                kernel_size=conv_kernels[0],
                stride=1,
            ),
            "conv_reluMixed": nn.ReLU(),
        }
        self.layers["convolution_mixed"] = {
            "conv3dM": nn.Conv3d(
                in_channels=hc[0],
                out_channels=hc[0]*2,
                kernel_size=conv_kernels[1],
                stride=1,
            ),
            "conv_reluM": nn.ReLU(),
        }
        self.layers["last_pool"] = {
            "last_pool": nn.AvgPool3d(
                kernel_size=last_kernel,
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
            "lineal_classifier": nn.Linear(latent_size**3, 1),
        }

        self.layers["decoder"] = {

        }
        self.layers["compressor"] = {
            "compressor1": nn.AvgPool3d(
                kernel_size=conv_kernels[0],
                stride=1,
            ),
            "compressor2": nn.AvgPool3d(
                 kernel_size=conv_kernels[1],
                 stride=1,
             ),
        }


        self.layers["default"] = {
            **self.layers["convolution_common"], # --> [3840, 10, 10, 10]
            **self.layers["convolution_mixed"],
            **self.layers["last_pool"],
            **self.layers["linear"], # --> [1, 216]
            **self.layers["classifier"],
        }

        self.criterions["default"] = ImgAndClassifierLoss()

    def forward(self, x, classify=False, reference=None):
        self.set_mode("default")
        x = self._forward(x, "convolution_common")

        x = self._forward(x, "convolution_mixed")
        x = self._forward(x, "last_pool")

        #print(x.shape)
        x = self._forward(x, "linear")
        if classify:
            c = self.classify(x, reference=reference)
        else:
            c = None
        x = x.reshape(1, self.data["latent_size"], self.data["latent_size"], self.data["latent_size"])
        if classify:
            return x, c
        return x

    def compress(self, x):
        with torch.no_grad():
            #print("compressing:", x.shape)
            x = self._forward(x, "compressor")
            x = self._forward(x, "last_pool")
            #print("compressed:", x.shape)
            return x

    def classify_norm(self, x, reference=None, return_diff=False):
        if reference is not None:
            x = x.reshape(1, self.data["latent_size"]**3)
            r = reference.to(x.dtype).reshape(1, self.data["latent_size"]**3).to(DEVICE)
            x = torch.nn.functional.normalize(x)
            r = torch.nn.functional.normalize(r)
            d = torch.sub(r, x)
        else:
            d = x

        c = self._forward(d, "classifier")
        c = c.reshape(1)
        if return_diff:
            d = d.reshape(1, self.data["latent_size"], self.data["latent_size"], self.data["latent_size"])
            return c, d
        return c

    def classify(self, x, reference=None, return_diff=False):
        if reference is not None:
            x = x.reshape(1, self.data["latent_size"]**3)
            r = reference.to(x.dtype).reshape(1, self.data["latent_size"]**3).to(DEVICE)

            max_x, min_x = torch.max(x), torch.min(x)
            max_r, min_r = torch.max(r), torch.min(r)
            #print("min/max")
            #print(max_x, min_x)
            #print(max_r, min_r)
            
            #print("Ranges")
            range_x = torch.abs(torch.sub(max_x, min_x))
            #print(range_x)
            range_r = torch.abs(torch.sub(max_r, min_r))
            #print(range_r)
            scaled_x = torch.multiply(torch.divide(torch.add(x, min_x), range_x),  range_r)

            d = torch.absolute(torch.sub(scaled_x, r))
        else:
            d = x

        c = self._forward(d, "classifier")
        c = c.reshape(1)
        if return_diff:
            d = d.reshape(1, self.data["latent_size"], self.data["latent_size"], self.data["latent_size"])
            return c, d
        return c
