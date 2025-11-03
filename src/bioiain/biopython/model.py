import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .chain import Chain



class Model(bp.Model.Model, BiopythonOverlayClass):
    child_class = Chain

