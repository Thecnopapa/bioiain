import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .atom import Atom

class Residue(bp.Residue.Residue, BiopythonOverlayClass):
    child_class = Atom