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
    def __init__(self):

        self.layers["default"] = {
            "linear1": nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-1]),
            "linear2": nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-1]),
            "linear3": nn.Linear(self.data["in_shape"][0], self.data["hidden_dims"][-1]),

        }








