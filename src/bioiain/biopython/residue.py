import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .atom import Atom, DAtom

class Residue(bp.Residue.Residue, BiopythonOverlayClass):
    child_class = Atom
    disordered_child_class = DAtom