import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .model import Model


class Structure(bp.Structure.Structure, BiopythonOverlayClass):
    child_class = Model
