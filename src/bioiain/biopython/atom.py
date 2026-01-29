import Bio.PDB as bp
from .base import BiopythonOverlayClass


class Atom(bp.Atom.Atom, BiopythonOverlayClass):
    child_class = None

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)


class BIAtom(BiopythonOverlayClass):
    child_class = None
    def __init__(self, data):
        if len(data) == 1:
            data = data[0]
        for k, v in data.items():
            if v == ".":
                data[k] = None
        self.type = data["group_PDB"]
        self.element = data["type_symbol"]
        self.atomnum = data["label_seq_id"]
        if self.atomnum is not None: self.atomnum = int(self.atomnum)
        self.name = data["label_atom_id"]
        self.resname = data["label_comp_id"]
        self.chain = data["auth_asym_id"]
        self.resnum = data["auth_seq_id"]
        if self.resnum is not None: self.resnum = int(self.resnum)
        self.id = (self.name,  self.resnum, self.chain)
        self.x = float(data["Cartn_x"])
        self.y = float(data["Cartn_y"])
        self.z = float(data["Cartn_z"])
        self.occupancy = float(data["occupancy"])
        self.b = float(data["B_iso_or_equiv"])
        self.coord = (self.x, self.y, self.z)
        self.disordered = False
        self.doppelgangers = []

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)





class DAtom(bp.Atom.DisorderedAtom, Atom):
    child_class = None

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for a in self.disordered_get_id_list():
            self[a] = Atom.cast(self.disordered_get(a))
        self.disordered_select(self.id)