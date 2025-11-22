import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .residue import Residue, DResidue

class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    child_class = Residue
    disordered_child_class = DResidue

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def _init(self, *args, **kwargs):
        pass