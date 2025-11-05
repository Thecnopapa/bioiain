import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .atom import Atom, DAtom

class Residue(bp.Residue.Residue, BiopythonOverlayClass):
    child_class = Atom
    disordered_child_class = DAtom

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)