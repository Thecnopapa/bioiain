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

    def _calculate_compactness(self, radius=10, plot=False, session=False):
        kdtree = self.ca_kdtree()
        print(kdtree)
        from src.bioiain.visualisation.plots import plasma

        if plot:
            from src.bioiain.visualisation.plots import fig3D, close, show, line
            fig, ax = fig3D()
            print(fig, ax)

        if session:
            from src.bioiain.visualisation.pymol import PymolScript
            script = PymolScript(name=f"compactness_{self.name()}", folder = self.folder())
            minimal = self.path(minimal=True)
            entity_name = script.load(minimal)
            print(script ,entity_name)


        max_frags = max([k["atom"].get_misc("fragment", 0) for k in kdtree])
        for n, k in enumerate(kdtree):
            if k["op"] != 1:
                continue
            # print(n, k)
            fragment = k["atom"].get_misc("fragment", 0)
            color = plasma(fragment, scale=max_frags)
            coord = k["coord"]
            atom = k["atom"]
            neighs = list(kdtree.radius(k["coord"], radius=radius)[0])
            # print(neighs)
            final_vector = np.array([0., 0., 0.])
            valid_nn = 0
            for nn in neighs:
                if n == nn:
                    # log("warning", "Same atom:", n, nn)
                    continue
                if (kdtree.atom_of(nn).get_misc("fragment", None) == fragment) and (kdtree.pos_of(nn) in [None, 1]):
                    # log("warning", "Same fragment:", fragment,  kdtree.atom_of(nn).get_misc("fragment", None),)
                    continue
                valid_nn += 1
                # ax.plot(*line(k["coord"], kdtree.coord_of(nn)), c=color)
                final_vector += np.array(vector(coord, kdtree.coord_of(nn)))



            if valid_nn > 0:
                final_vector /= valid_nn
                final_vector *= -1
                compactness = length(final_vector)
                # print(final_vector, compactness)
                ccol = plasma(compactness, scale=10)
                hexccol = plasma(compactness, scale=10, as_pymol_hex=True)
                vector_end = coord + final_vector
                if plot:
                    ax.scatter(*coord, c=ccol, s=valid_nn + 1)
                    ax.plot(*line(coord, vector_end), c=ccol)
                if session:
                    r1 = f"({entity_name} and i. {atom.resnum} and c. {atom.complex})"
                    #print(r1)
                    a1 = f"({r1} and n. ca)"
                    script.line(name="compactness", sele1=a1, coord2=vector_end)
                    script.color(r1, color=hexccol)


        if plot:
            show()
            close(fig)

        if session:
            script.compile()
            script.execute()













