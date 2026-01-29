import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .residue import Residue, DResidue
from ..utilities.sequences import d3to1
from ..utilities.logging import log


class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    child_class = Residue
    disordered_child_class = DResidue

    def _init(self, *args, **kwargs):
        self.data["sequence"] = None


    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)


    def get_sequence(self, force=False):
        from .imports import read_mmcif

        atoms = read_mmcif(self.paths["self"], subset=["_atom_site"])("_atom_site")
        atoms = [a for a in atoms.values() if a["label_atom_id"] == "CA"]

        if self.data["sequence"] is not None and not force:
            return self.data["sequence"]
        self.data["sequence"] = ""
        for atom in atoms:
            try:
                #print(atom["label_comp_id"])
                r = d3to1[atom["label_comp_id"]]
                #print(r)
                self.data["sequence"] += r
            except KeyError as e:
                log("debug", f"Residue {e} not recognised")
                #input("press enter")
        return self.data["sequence"]
