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
    def __init__(self, *args, optimiser, use_original=True, range=4, **kwargs):
        self.o_lrs = [p["lr"] for p in optimiser.param_groups]
        self.use_original = use_original
        self.loss_list = []
        super().__init__(optimiser, *args, **kwargs)

    def get_lr(self):
        print("LRS: getting lrs")
        print(self.lrs)
        return torch.Tensor(np.array(self.lrs))


    def step(self, running_loss=None):
        print("LRS: stepping...")
        print("Current loss:", running_loss)

        if running_loss is None:
            return self.o_lrs

        if len(self.loss_list) == 0:
            self.loss_list.append(running_loss)
        print("Previous loss", self.loss_list[-1])
        self.lrs = []
        for p, olr in zip(self.optimizer.param_groups, self.o_lrs):
            if self.use_original:
                old_log = math.log(olr, 10)
                new_log = old_log - ((1-running_loss)*4) + 2
            else:
                old_log = math.log(p["lr"], 10)
                delta = self.loss_list[-1] / running_loss
                new_log = old_log - delta + 1 # -0.5


            print("OLD_LOG", old_log)

            print("NEW_LOG", new_log)
            new_lr = 10 ** new_log
            print("NEW_LR", new_lr)
            self.lrs.append(new_lr)
            p["lr"] = torch.Tensor(np.array([new_lr]))

        self.loss_list.append(running_loss)
        print(self.lrs)
        return self.lrs




class CustomLoss(object):
    pass




class VQLoss(CustomLoss):
    def __init__(self, weights=None):
        self.CEL = nn.MSELoss()

    def encoder_loss_vq_vae(self, o, tensor, commitment=0.25):
        loss = self.CEL(o, tensor)
        loss = loss * (1 + commitment)
        return loss

    def encoder_loss(self, o, tensor, origin, commitment=0.25):
        t_loss = self.CEL(o, tensor)
        o_loss = self.CEL(o, tensor)

        loss = o_loss-t_loss

        loss = loss * (1 + commitment)
        return loss

    def decoder_loss(self, o, tensor):
        loss = self.CEL(o, tensor)
        return loss

    def __call__(self, e_loss, d_loss):

        loss = e_loss + d_loss

        return loss




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

