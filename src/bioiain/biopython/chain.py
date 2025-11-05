import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .residue import Residue

class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    child_class = Residue

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def init(self, *args, **kwargs):
        pass