import os, json

from src.bioiain.utilities import *
from src.bioiain.utilities.sequences import d3
from src.bioiain.utilities.exceptions import *
from src.bioiain.machine.embeddings import *

from torch import Tensor


device = "cpu"





class ALEPHEmbedding(ResidueEmbedding):
    param_names = ["len_i", "len_j", "len_ij", "angle_ij", "dihedral_ij", "theta_i", "theta_j"]

    def __init__(self, *args, cvector=None, modulo_norm=2.4, max_dist=20, **kwargs):
        self.cvector = cvector
        self.modulo_norm = modulo_norm
        self.max_dist = max_dist
        super().__init__(name=self.cvector.full_id(), residue=self.cvector.res2)


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
    param_names = ["len_i", "len_j", "len_ij", "angle_ij", "dihedral_ij", "theta_i", "theta_j", "bfactor", "contactability", "dist_to_ligand"]
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def _generate(self) -> list:
        cv = self.cvector
        e = super()._generate()

        # print(cv.chain, cv.closest.chain, cv.chain != cv.closest.chain, cv.closest_opn, cv.closest_opn not in [None, 1])
        if cv.chain != cv.closest.chain or cv.closest_opn not in [None, 1]:
            is_contact = 1
        else:
            is_contact = 0
        b = min(cv.res2.bfactor(), 100) /100
        dl = cv.dist_to_lig
        dl = min(1, dl / self.max_dist)


        e.extend([b, is_contact, dl])

        return e



class ALEPHProteinEmbedding(ProteinEmbedding):
    residue_embedding_class = ALEPHEmbedding

    def check_aleph(self, *args, vc_mode=None, in_place=True, **kwargs) -> bool:

        if not self.entity.has_flag("no_atoms", True):
            try:
                self.entity = self.entity.fragment(in_place=in_place)
            except ALEPHError:
                raise

            if self.entity.has_flag("missing_side_chains"):
                raise NoEmbeddingForThisProtein()

            if self.entity.data["fragments"]["n_fragments"] <= 1:
                raise NoEmbeddingForThisProtein()

            cvectors =  self.entity.cvectors(vc_mode=vc_mode)
            cvmatrix =  self.entity.cvmatrix(vc_mode=vc_mode)  # Not used but calculates closest neighbours
            if cvmatrix is None:
                raise NoEmbeddingForThisProtein

    def __init__(self, *args,  **kwargs):
        super().__init__(*args, **kwargs)

    def generate(self, *args, **kwargs) -> Tensor:
        self.check_aleph(*args, **kwargs)
        return super().generate(*args, **kwargs)

    def _generate(self, *args, **kwargs) -> list:
        assert self.residue_embedding_class is not None
        e = []
        seq = ""
        for n, cv in enumerate(self.entity.cvectors()):
            try:
                e.append(self.residue_embedding_class(*args, cvector=cv, **kwargs).tensor())
                seq += d3(cv.res2.resname)[0]
            except NoEmbeddingForThisResidue:
                seq += "-"
                self.missing_indexes.append(n)
        self.sequence = seq
        return e
