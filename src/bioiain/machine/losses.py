import os, json
import numpy as np

from ..utilities.logging import *
from ..utilities.parallel import *

import torch
import torch.nn as nn


from .datasets import Item

from . import DEVICE

from ..utilities.exceptions import *



class CustomLoss(object):
    pass





class CustomWeighted(CustomLoss):
    def __init__(self, weights=None, addto1=False):
        log(1, "Using label weights:")

        if weights is not None:

            if type(weights) is dict:
                weights = weights.values()
            w = np.array(list(weights))
            print(w)
            if addto1:
                w = w / w.sum()
            else:
                w = w / w.max()

            self.weight = torch.Tensor(w)

        self.CEL = nn.MSELoss()
        log(2, self.weight)

    def __call__(self, o, item):
        true_index = item.li
        true_tensor = item.lt

        loss = self.CEL(o, true_tensor)
        
        pred_index = torch.max(o, dim=0)[1]

        if self.weight is not None:
            pred_weight = self.weight[pred_index]
            loss *= pred_weight

        return loss




class CustomHalfHalf(CustomLoss):
    def __init__(self, weights=None):
        log(1, "Using label weights:")

        if type(weights) is dict:
            weights = weights.values()
        w = np.array(list(weights))
        w = w
        w = w / w.sum()
        w = torch.Tensor(w)
        log(2, w)
        self.CEL = nn.MSELoss()
        self.weight = w

    def __call__(self, o, item):
        true_index = item.li
        true_tensor = item.lt
        pred = torch.max(o, dim=0)[1]
        if len(item.lt) == 10:
            if item.li < 5:
                true_tensor[:5] = 0.5
            else:
                true_tensor[5:] = 0.5
            true_tensor[item.li:item.li+1] = 1.

        #print(true_tensor)
        weighted_out = o# * self.weight
        #print("Weighted out:", weighted_out)
        loss = self.CEL(o, true_tensor)
        #if true_index < 5 == torch.max(o, dim=0)[1] < 5:
        #    loss *= 0.5
        loss *= self.weight[pred]
        return loss

