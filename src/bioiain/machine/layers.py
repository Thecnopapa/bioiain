import os, json
import numpy as np

from ..utilities.logging import *
from ..utilities.parallel import *

import torch
import torch.nn as nn


from .datasets import Item

from . import DEVICE

from ..utilities.exceptions import *
from ..utilities.maths import *




class Codebook(nn.Module):
    def __init__(self, n_tokens=20, latent_dims=2, commitment=0.25):
        super().__init__()
        log(1, "Initialising codebook...")
        self.commitment = commitment
        self.n_tokens = n_tokens
        self.latent_dims = latent_dims
        self.codebook = nn.Embedding(n_tokens, latent_dims)
        self.codebook.weight.data.uniform_(-1/self.n_tokens, 1/self.n_tokens)
        self.MSE = nn.MSELoss()
        self.last_loss = None
        self.last_index = None

        log(2,"Codebook:", self.codebook)
        #print("weight:")
        #print(self.codebook.weight)



    def forward(self, x):
        #print("forward")
        #print(x)
        # Calculate distances between z and the codebook embeddings |a-b|²
        distances = (
            torch.sum(x ** 2, dim=-1, keepdim=True)                 # a²
            + torch.sum(self.codebook.weight.t() ** 2, dim=0, keepdim=True)  # b²
            - 2 * torch.matmul(x, self.codebook.weight.t())        # -2ab
        )
        #print(distances)
        closest_index = torch.argmin(distances, dim=1)
        closest_tensor = self.codebook(closest_index)

        #print("closest:", closest_index, "distance:", distances[0][closest_index])
        #print(closest_tensor)
        loss = self.MSE(closest_tensor, x.detach()) + self.commitment * self.MSE(closest_tensor.detach(), x)

        closest_tensor = x + (closest_tensor - x).detach()
        self.last_loss = loss
        #print(closest_index)
        self.last_index = closest_index

        return closest_tensor












