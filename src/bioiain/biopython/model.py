import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .chain import Chain



class Model(bp.Model.Model, BiopythonOverlayClass):
    child_class = Chain

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

