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
