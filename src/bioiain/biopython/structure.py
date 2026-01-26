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
        from src.bioiain.symmetries import MissingCrystalError
        try:
            from ..symmetries import parse_crystal_card, calculate_parameters
            self.data["crystal"] =  parse_crystal_card(self.paths["original"])
            self.data["params"] = calculate_parameters(self.data["crystal"])
        except MissingCrystalError:
            self.data["crystal"] = {}
            self.data["crystals"] = []
            self.data["params"] = None
            return None
        self.pass_down()
        self.export()
        self.data["crystals"] = self.get_crystals()
        self.export()
        return self.data["crystals"]



    def init_all(self):
        self.init_crystal()
        self.pass_down()
        self.export()

    def get_crystals(self, first_model_only=True):
        from ..symmetries import Crystal
        if first_model_only:
            c = Crystal.cast(self.get_list()[0].copy(), export=True)
            self.paths["crystal_folder"] = c.paths["export_folder"]
            self.data["crystals"] = [c.get_name()]
        else:
            models = self.get_list()
            cs = [Crystal.cast(m.copy(), export=True) for m in models]
            self.paths["crystal_folder"] = cs[0].paths["export_folder"]
            self.data["crystals"] = [c.get_name() for c in cs]
        self.export()
        return self.data["crystals"]
