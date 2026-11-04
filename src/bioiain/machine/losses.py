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




class CustomLRS(torch.optim.lr_scheduler.LRScheduler):
    def __init__(self, optimiser, min=-10, max=-1, **kwargs):
        self.o_lrs = [p["lr"] for p in optimiser.param_groups]
        self.loss_list = []
        self.min = min
        self.max = max
        self.kwargs = kwargs
        super().__init__(optimiser)

    def __repr__(self):
        return f"<{self.__class__.__name__}: {self.optimizer} (x{self.factor})>"

    def get_lr(self):
        log(2, "Getting lrs:", [f"{l:.3e}" for l in self.lrs])
        return torch.Tensor(np.array(self.lrs))

    def current(self):
        try:
            return self.lrs
        except AttributeError:
            return self.o_lrs

    def last(self):
        try:
            return self._last_lrs
        except AttributeError:
            return self.o_lrs

    def _step(self, running_loss, **kwargs):
        raise NotImplementedError()

    def step(self, running_loss):
        log(2, f"Stepping: {self}")
        if running_loss is None:
            return self.o_lrs
        log(3, f"Current loss: {running_loss:.3e}")



        if len(self.loss_list) == 0:
            self.loss_list.append(running_loss)
        log(3, f"Previous loss {self.loss_list[-1]:.3e}")
        if hasattr(self, "lrs"):
            self._last_lrs = self.lrs
        self.lrs = []

        self._step(running_loss, **self.kwargs)

        self.loss_list.append(running_loss)
        log(2, f"Updating LRs:", [f"{l:.3e}" for l in self.lrs])
        super().step()
        return self


    def clamp(self, val):
        return max(min(self.max, val), self.min)



class PlateauLRS(torch.optim.lr_scheduler.ReduceLROnPlateau):
    def current(self):
        return self.get_last_lr()

    def last(self):
        return self.get_last_lr()


class EpochLRS(CustomLRS):
    def _step(self, running_loss, frequency=1, delta=-0.001, **kwargs):
        if len(self.loss_list) % frequency != 0:
            log(4, f"Not stepping in this epoch")
            return
        for n, (p, olr) in enumerate(zip(self.optimizer.param_groups, self.last())):
            old_log = math.log(olr, 10)
            log(4, f"OLD_LOG {old_log:.3f}")

            new_log = old_log + delta
            log(4, f"NEW_LOG: {new_log:.3f}")

            new_log = self.clamp(new_log)
            log(4, f"CAPPED NEW_LOG: {new_log:.3f}")

            new_lr = 10 ** new_log
            log(4, f"NEW_LR {new_lr:.3e}")

            self.lrs.append(new_lr)
            p["lr"] = float(new_lr)

class DeltaLRS(CustomLRS):
    def _step(self, running_loss, factor=1, threshold=0.5, smoothing=5, freedom=0.5):
        for n, (p, olr) in enumerate(zip(self.optimizer.param_groups, self.last())):
            old_log = math.log(olr, 10)

            smoothed_loss = 0
            smoothed_n = 0
            for i in range(smoothing):
                try:
                    smoothed_loss += self.loss_list[-i-1]
                    smoothed_n += 1
                except IndexError:
                    break
            smoothed_loss /= smoothed_n
            log(4, f"Smoothed loss: {smoothed_loss:.3e}, {smoothed_n}")
            ratio = smoothed_loss / running_loss
            log(4, f"RATIO: {ratio:.3f}")


            if running_loss < threshold:
                log(4, "Below threshold!")
                if ratio >= 1:
                    delta = -(ratio - 1)*factor
                else:
                    delta = 0


            elif ratio >= 2:
                #delta = delta = -(ratio - 1) /2
                delta = 0
            elif ratio >= 1:
                delta = 0
            elif ratio <= 0.5:
                delta = -(0.5 - ratio)*factor
            else:
                delta = (1 - ratio)*factor

            # + DELTA -> + LR

            log(4, f"DELTA: {delta:.3f}")

            delta = max(min(delta, freedom), -freedom)

            log(4, f"CAPPED DELTA: {delta:.3f}")


            log(4, f"OLD_LOG {old_log:.3f}")
            new_log = old_log + delta
            log(4, f"NEW_LOG: {new_log:.3f}")

            new_log = self.clamp(new_log)
            log(4, f"CAPPED NEW_LOG: {new_log:.3f}")

            new_lr = 10 ** new_log
            log(4, f"NEW_LR {new_lr:.3e}")

            self.lrs.append(new_lr)
            p["lr"] = float(new_lr)





class CustomLoss(object):
    pass




class VQLoss(CustomLoss):
    def __init__(self, weights=None):
        self.MSE = nn.MSELoss()
        self.CSL = nn.CosineEmbeddingLoss(margin=0.)

    def encoder_loss_vq_vae(self, o, tensor, commitment=0.25):
        loss = self.MSE(o, tensor)
        loss = loss * (1 + commitment)
        return loss

    def encoder_loss(self, l, token, origin, commitment=0.25):
        token = token.to(DEVICE)
        origin = origin.to(DEVICE)
        #print()
        #print("L", l)
        #print("T", token)
        #print("O", origin)
        t_l = torch_distance(l, token)
        o_t = torch_distance(token, origin)
        o_l = torch_distance(l, origin)
        y = torch.Tensor([-1])[0]
        z = torch.zeros(1).to(DEVICE)

        #print("TL", t_l)
        #print("OT", o_t)
        #print("OL", o_l)
        #print("Y", y)

        #loss = self.CSL(t_l, o_l, y)
        #print(loss)
        #print("LOSS", l_t.item(), "-", (o_l-o_t).item())


        #loss = t_l-(o_l-o_t)
        loss = (t_l-o_l+o_t) / t_l




        loss = torch.max(z, loss * (1 + commitment))
        #print(loss)
        return loss

    def decoder_loss(self, o, tensor):
        loss = self.MSE(o, tensor)
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








class ImgAndClassifierLoss(object):
    def __init__(self):
        self.MSE = nn.MSELoss()

    def __call__(self, out_i, img, out_c=None, c=None, return_all=False):
        #print(out_i, img, out_c, c)
        #print("OUT_I:", out_i)
        #print("IMG:", img)
        i_loss = self.MSE(out_i, img)

        c_loss = None
        try:

            assert (out_c is not None) and (c is not None), "No class provided"
            print("\nI LOSS:", i_loss)

            c_loss = self.MSE(out_c, c)
            print("C LOSS:", c_loss)
            # if c_loss < 0.001:
            #     raise LossIsZero()
            # multiplier = torch.divide(torch.tensor([1]).to(DEVICE), c_loss)
            multiplier = c_loss

            print("MULTIPLIER:", multiplier)
            loss = torch.multiply(i_loss, multiplier)
            print("NEW LOSS", loss)

        except (AssertionError, LossIsZero):
            loss = i_loss

        #print(str(loss.detach().cpu().numpy()))
        if str(loss.detach().cpu().numpy()) == "nan":
            raise Exception("Loss is None")
        if return_all:
            return loss, i_loss, c_loss
        return loss



class TwoCubeLoss(object):
    def __init__(self):
        self.MSE = nn.MSELoss()

    def __call__(self, cube1, cube1label, cube2=None, cube2label=None, return_all=False):


        #print(cube1.shape)
        #print(cube1label.shape)
        #print(cube2.shape)
        #print(cube2label.shape)
        loss1 = self.MSE(cube1, cube1label)
        loss2 = self.MSE(cube2, cube2label)

        loss = torch.add(loss1, loss2)

        if loss.isnan():
            print(loss, loss1, loss2)
            if loss1.isnan():
                print(cube1)
                print(cube1label)
            if loss2.isnan():
                print(cube2)
                print(cube2label)
            raise Exception(f"Loss is {loss} ({type(loss)})")


        if return_all:
            return loss, loss1, loss2

        return loss
