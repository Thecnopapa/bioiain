import os, sys, json
import numpy as np

from ..utilities import *
from ..utilities.exceptions import *
from ..utilities.maths import *


from ..utilities.space_groups import dictio_space_groups


class SpaceGroup(object):
    def __init__(self, sp, key_is_str=False):

        group = None
        if (type(sp) is int) or key_is_str:
            n = sp
            group = dictio_space_groups[sp]
        else:
            for n, group in dictio_space_groups.items():
                if sp.strip() in [group.get("symbol", None), group.get("xHM_symbol", None), group.get("short_symbol", None), group.get("hall_symbol", None)]:
                    break

        if group is None:
            raise SpaceGroupNotFound(sp)
        self.group = group
        self.n = n
        if type(n) is str:
            log("warning", f"Space group does not match an International Table entry: {n}")
        self.name = self.group["xHM_symbol"]












class UnitCell(object):
    def __init__(self, a, b, c, alpha, beta, gamma, space_group=None, z=None):
        self.a = a
        self.b = b
        self.c = c
        self.alpha = alpha
        self.beta = beta
        self.gamma = gamma
        self.z = None
        self.space_group = space_group

    def __repr__(self):
        pstr = " ".join([f"{p:.2f}" for p in self.params])
        return f"<bi.{self.__class__.__name__}[{pstr}]"


    def params(self):
        return [self.a, self.b, self.c, self.alpha, self.beta, self.gamma]

    def reduce(self, mode="niggli"):
        pass

    def is_reduced(self, mode="niggli"):
        pass