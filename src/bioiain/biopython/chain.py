import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .residue import Residue

class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    child_class = Residue