import os, json
import numpy as np

from ..utilities.logging import *
from ..utilities.parallel import *

import torch
import torch.nn as nn


from .datasets import Item

from . import DEVICE

from ..utilities.exceptions import *




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

