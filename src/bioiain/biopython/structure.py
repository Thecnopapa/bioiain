import Bio.PDB as bp
from .base import BiopythonOverlayClass


class Structure(bp.Structure.Structure, BiopythonOverlayClass):
    pass