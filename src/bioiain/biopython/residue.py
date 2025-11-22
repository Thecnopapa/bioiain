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

class DResidue(bp.Residue.DisorderedResidue, Residue):
    child_class = Atom
    disordered_child_class = DAtom

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for a in self.disordered_get_id_list():
            self[a] = Residue.cast(self.disordered_get(a))
        self.disordered_select(self.id)