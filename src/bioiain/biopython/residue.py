import Bio.PDB as bp
from .base import BiopythonOverlayClass


class Residue(bp.Residue.Residue, BiopythonOverlayClass):
    pass