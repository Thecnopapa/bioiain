import Bio.PDB as bp
from .base import BiopythonOverlayClass


class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    pass