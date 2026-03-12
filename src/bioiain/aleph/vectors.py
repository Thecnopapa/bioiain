import os, json
import numpy as np

from ..utilities.logging import log
from ..utilities.exceptions import *








class CVector(object):
    def __init__(self, res1=None, res2=None, res3=None):
        self.res1 = res1
        self.res2 = res2
        self.res3 = res3
        self.residues = (res1, res2, res3)
        self.chain = res2.chain

        self.length = None
        self.start = None
        self.end = None

        self.calculate()


    def __repr__(self):
        return f"<bi.CVector c.{self.chain} (>{self.res1.resnum}-{self.res2.resnum}-{self.res3.resnum}>)>"


    def _centroid(self, *coords):
        coords = np.array(coords)
        coords = coords.sum(axis=0)
        coords = coords/len(coords)
        return coords



    def calculate(self):

        self.start = self._centroid([r.ca.coord for r in self.residues])
        self.end = self._centroid([r.o.coord for r in self.residues])

        self.length = self.end - self.start
        return self




