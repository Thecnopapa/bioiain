import os, json

from ..utilities.exceptions import *
from .entity import BIEntity
from .residue import BIResidue

class BIChain(BIEntity):
    child_class = BIResidue
    extension = "chain"
    level = "chain"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.paths["sub_folder"] = "chains"

    def find_id(self, method="first"):

        if method == "first":
            chain_id = self.atoms()[0].chain
        else:
            chain_id = method
        return chain_id


    def id(self):
        return self.data["info"]["chain_id"]



    def set_chain_id(self, chain_id=None, complex=False):
        if chain_id is None:
            chain_id = self.find_id()
        chain_id = str(chain_id)
        if len(chain_id) != 1:
            chain_id = chain_id.replace("_", "-")

            log(f"warning", "CHAIN ID is longer than 1: {chain_id}")
            self.set_flag("unconventional_chain_id", True)
        self.data["info"]["chain_id"] = chain_id
        if self.has_flag("has_chain_id", True):
            old_name = self.name().split("_")
            for n, o in enumerate(old_name):
                if o == self.id():
                    old_name[n] = self.id()
            new_name = "_".join(old_name)
            self.set_name(new_name, append=False)
        else:
            self.set_name(chain_id, append=True)
            self.set_flag("has_chain_id", True)

        for a in self.all_atoms():
            a.chain = chain_id
            if complex:
                a.complex = chain_id
        return self.id()



    @classmethod
    def from_atoms(cls, atoms, code=None, chain_id=None, overwrite_complex=False, **kwargs):
        self = super().from_atoms(atoms, code, **kwargs)
        self._atoms = atoms
        self.set_chain_id(chain_id, complex=overwrite_complex)
        self.sequence()
        return self

    @classmethod
    def from_file(cls,*args, chain_id=None, overwrite_complex=False, **kwargs):
        self = super().from_file(*args, **kwargs)
        self._atoms = self.atoms(chain=chain_id, hetatm=True, disordered=True)
        self.set_chain_id(chain_id, complex=overwrite_complex)
        self.sequence()
        return self

