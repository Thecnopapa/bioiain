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




def compactness_mpl(kdtree, radius=10):
    from src.bioiain.visualisation.plots import fig3D, close, show, plasma, line

    fig, ax = fig3D()
    print(fig, ax)

    max_frags = max([k["atom"].get_misc("fragment", 0) for k in kdtree])
    for n, k in enumerate(kdtree):
        if k["op"] != 1:
            continue
        print(n, k)
        fragment = k["atom"].get_misc("fragment", 0)
        color = plasma(fragment, scale = max_frags)
        coord = k["coord"]
        neighs = list(kdtree.radius(k["coord"], radius=radius)[0])
        print(neighs)
        final_vector = np.array([0., 0., 0.])
        valid_nn = 0
        for nn in neighs:
            if n == nn:
                log("warning", "Same atom:", n, nn)
                continue
            if kdtree.atom_of(nn).get_misc("fragment", None) == fragment:
                log("warning", "Same fragment:", fragment,  kdtree.atom_of(nn).get_misc("fragment", None),)
                continue
            valid_nn += 1
            #ax.plot(*line(k["coord"], kdtree.coord_of(nn)), c=color)
            final_vector += np.array(vector(coord, kdtree.coord_of(nn)))
        if valid_nn > 0:
            final_vector /= valid_nn
            final_vector *= -1
            compactness = length(final_vector)
            print(final_vector, compactness)
            ccol = plasma(compactness, scale = 10)
            ax.scatter(*coord, c=ccol)
            ax.plot(*line(coord, coord+final_vector), c=ccol)










    show()
    close(fig)






class CompactStructure(FragmentedStructure):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def ca_kdtree(self, force=False, **kwargs):
        if force or getattr(self, "_kdtrees", {}).get("ca", None) is None:
            KDT(self, mode="ca", auto_parse_symmetry=True, **kwargs)
        return self._kdtrees["ca"]


    def _calculate_compactness(self, radius=10):

        kdtree = self.ca_kdtree()
        print(kdtree)
        # for n, k in enumerate(kdtree):
        #     print(n, k)
        #     print(kdtree.radius(k["coord"], radius=radius))

        compactness_mpl(kdtree)











