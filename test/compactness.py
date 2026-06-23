import os, sys, json

sys.path.append('..')

from src.bioiain.base import BIEntity
from src.bioiain.aleph import FragmentedStructure
import polars as pl
import numpy as np
from src.bioiain.utilities.exceptions import *
from src.bioiain import log
from src.bioiain.utilities import *
from src.bioiain.utilities.maths import *
from src.bioiain.utilities.kdtree import KDT

class CompactStructure(FragmentedStructure):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def ca_kdtree(self, force=False, **kwargs):
        if force or getattr(self, "_kdtrees", {}).get("ca", None) is None:
            KDT(self, mode="ca", auto_parse_symmetry=True, **kwargs)
        return self._kdtrees["ca"]


    def _calculate_compactness(self):

        kdtree = self.ca_kdtree()
        print(kdtree)
        for k in kdtree:
            print(k)











