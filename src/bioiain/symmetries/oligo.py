from .space_groups import dictio_space_groups
import Bio.PDB as bp

from .operations import *
from ..utilities.logging import log
from ..biopython import Chain, Model





class Crystal(Model):
    def init(self, *args, **kwargs):
        self.data["min_monomer_length"] = None
        self.data["oligomer_levels"] = None

        self.monomers = None
        self.ligands = None

    def set_params(self,params, min_monomer_length, oligomer_levels):
        self.data["params"] = params
        self.data["min_monomer_length"] = min_monomer_length
        self.data["oligomer_levels"] = oligomer_levels

    def process(self):
        self._identyfy_main_elements()
        self._cast_main_elements()

    def _identyfy_main_elements(self):
        if self.data["min_monomer_length"] is None:
            log("error", "Crystal: missing param: min_monomer_length", raise_exception=True)
        self.monomers = []
        self.ligands = []
        for chain in self.get_chains():
            c_len = len(chain)
            print(chain, c_len)
            if c_len <= self.data["min_monomer_length"]:
                self.ligands.append(chain)
            else:
                self.monomers.append(chain)


    def _cast_main_elements(self):
        for n, mon in enumerate(self.monomers):
            m = Monomer.cast(mon)
            self.monomers[n] = m
        for n, lig in enumerate(self.ligands):
            l = Ligand.cast(lig)
            self.ligands[n] = l
        log("debug", "Monomers: {}, Ligands: {}".format(self.monomers, self.ligands))


    def _find_oligomers(self):

        for n in self.data["oligomer_levels"]:
            log("debug" "Oligomer level: {}".format(n))




    def _regenerate_crystal(self, chain):
        pass



class CrystalElement(Chain):
    pass


class Monomer(CrystalElement):
    pass



class Ligand(CrystalElement):
    pass














