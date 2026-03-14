import os, json
import numpy as np

from . import CVMatrix
from ..base import BIStructure, BIChain, BIResidue
from ..utilities import *





class Fragment(BIChain):
    child_class = BIResidue
    extension = "fragment"

    def __init__(self, *args, fragment_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.paths["sub_folder"] = "fragments"


    @classmethod
    def from_atoms(cls, atoms, code=None, chain_id=None, fragment_id=None, **kwargs):
        self = super().from_atoms(atoms, code, chain_id, **kwargs)
        self.data["info"]["fragment_id"] = f"F{fragment_id}"
        self.set_name(fragment_id, append=True)
        for a in self.atoms():
            a.set_misc("fragment", fragment_id)
        return self






class FragmentedStructure(BIStructure):
    child_class = Fragment
    extension = "fstructure"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.paths["sub_folder"] = "fragments"
        self._cvectors = None



    def fragment_with_aleph(self):
        from .core.ALEPH import annotate_pdb_model_with_aleph
        self.export()
        print("### ALEPH start ###")
        graph, _, _, _, _ = annotate_pdb_model_with_aleph(
            self.paths["self"],
            weight="distance_avg",
            strictness_ah=0.45,
            strictness_bs=0.20,
            peptide_length=3,
            write_pdb=False,
        )
        print("### ALEPH end ###")

        log(1, "Fragmented:", self)



        fragments = []
        fragmented_ids = graph.vs

        for n, fraglist in enumerate(fragmented_ids):
            cvs = fraglist["reslist"]
            reslist = self.residues()
            target_res= []
            for fres in cvs:
                _, model, chain, full_id, resname = fres

                for res in reslist.copy():
                    if res.resname == resname and res.chain == chain and res.resnum == full_id[1]:
                        #print(res, fres)
                        target_res.append(res)
                        reslist.remove(res)
                        break
            log(2, f"fragment {n}: {len(target_res)}")
            fatoms = []
            for res in target_res:
                fatoms.extend(res.atoms)
            chain = list(set([r.chain for r in fatoms]))
            assert len(chain) == 1
            chain = chain[0]
            fragment = Fragment.from_atoms(fatoms, code=self.code(), chain_id=chain, fragment_id=n, parent=self)
            fragment.export()
            fragments.append(fragment)

        self._fragments = fragments
        return fragments


    def _calculate_cvectors(self):
        print("CALCULATING cvectors")
        from ..aleph.vectors import CVector
        residues = self.residues()
        n_res = len(residues)
        cvector_list = []
        for n, res in enumerate(residues):
            if n == 0 or n == n_res -1:
                continue

            cvector = CVector(residues[n-1], res, residues[n+1])
            cvector_list.append(cvector)

        self._cvectors = cvector_list
        return cvector_list

    def map_cvectors(self):
        if self._cvectors is None:
            self._cvectors = self._calculate_cvectors()
        matrix = CVMatrix(self._cvectors)
        return matrix








