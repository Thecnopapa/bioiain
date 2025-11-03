import Bio.PDB as bp
from .base import BiopythonOverlayClass


class Model(bp.Model.Model, BiopythonOverlayClass):
    pass
