import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .model import Model


class Structure(bp.Structure.Structure, BiopythonOverlayClass):
    child_class = Model

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def init_crystal(self):
        from ..symmetries import parse_crystal_card, calculate_parameters
        self.data["crystal"] =  parse_crystal_card(self.paths["original"])
        self.data["params"] = calculate_parameters(self.data["crystal"])


    def init_all(self):
        self.init_crystal()