import os, json

from src.bioiain.utilities import *
from src.bioiain.utilities.parallel import avail_cpus
from src.bioiain.utilities.sequences import d3
from src.bioiain.utilities.exceptions import *
from src.bioiain.machine.embeddings import *

import torch
from torch import Tensor


device = "cpu"





class ALEPHEmbedding(ResidueEmbedding):
    def __init__(self, cvector, modulo_norm=2.4, max_dist=20):
        self.cvector = cvector
        self.modulo_norm = modulo_norm
        self.max_dist = max_dist
        super().__init__(name=self.cvector.full_id(), residue=self.cvector.res2)
        self.param_names.extend(["len_i", "len_j", "len_ij", "angle_ij", "dihedral_ij", "theta_i", "theta_j"])


    def _generate(self) -> list:
        cv = self.cvector
        i = cv
        j = cv.closest
        i_j = cv.closest_vp

        len_i = min(1, i.d / self.modulo_norm)
        len_j = min(1, j.d / self.modulo_norm)
        len_i_j = min(1, i_j.d / self.max_dist)
        angle_i_j = min(1, i_j.a / 360)
        da = min(1, i_j.da / 360)
        t1 = min(1, i_j.t1 / 360)
        t2 = min(1, i_j.t2 / 360)

        e = [len_i, len_j, len_i_j, angle_i_j, da, t1, t2]
        return e


class ExpandedALEPHEmbedding0(ALEPHEmbedding):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.param_names.extend(["contactability", "dist_to_ligand"])

    def _generate(self) -> list:
        cv = self.cvector
        e = super()._generate()

        # print(cv.chain, cv.closest.chain, cv.chain != cv.closest.chain, cv.closest_opn, cv.closest_opn not in [None, 1])
        if cv.chain != cv.closest.chain or cv.closest_opn not in [None, 1]:
            is_contact = 1
        else:
            is_contact = 0

        dl = cv.dist_to_lig
        dl = min(1, dl / self.max_dist)


        e.extend([is_contact, dl])

        return e






class CVEmbedding():
    """
    Default CV embedding: 4 params [l(i), l(j) , a(ij), d(ij)].
    Normalised to 1 (2.4A, 10A by default).
    Virtual center is CV start.
    """
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)


    def _cvectors_to_embedding(self, cvectors, modulo_norm, max_dist, **kwargs):
        e = []
        seq = ""
        for cv in cvectors:
            i = cv
            j = cv.closest
            i_j = cv.closest_vp

            len_i = min(1, i.d / modulo_norm)
            len_j = min(1, j.d / modulo_norm)
            len_i_j = min(1, i_j.d/max_dist)
            angle_i_j = min(1, i_j.a / 360)
            rn, ri = d3(cv.resname)
            seq += rn

            e.append([len_i, len_j, angle_i_j, len_i_j])
        return e, seq

    def generate_embedding(self, *args, modulo_norm=2.4, max_dist=10, vc_mode=None, **kwargs):

        try:
            self.entity = self.entity.fragment()
            frag = self.entity
        except ALEPHError:
            return None

        if self.entity.has_flag("missing_side_chains"):
            return None

        if frag.data["fragments"]["n_fragments"] <= 1:
            return None
        cvectors = frag.cvectors(vc_mode=vc_mode)
        cvmatrix = frag.cvmatrix(vc_mode=vc_mode) # Not used but calculates closest neighbours
        if cvmatrix is None:
            return None

        e, seq = self._cvectors_to_embedding(cvectors,modulo_norm=modulo_norm, max_dist=max_dist, **kwargs)

        final_e = []
        final_seq = []
        for emb, seq in zip(e, seq):
            if any([ee is None for ee in emb]):
                continue
            final_e.append(emb)
            if type(seq) is list:
                seq = "".join(seq)
            final_seq.append(seq)

        final_e = torch.Tensor(final_e)
        torch.save(final_e, self.path)
        #print(e, e.shape, len(seq))
        self.sequence = final_seq
        self.length = len(self.sequence)
        self.exists = True
        return self.path




