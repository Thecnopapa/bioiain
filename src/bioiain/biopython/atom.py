import Bio.PDB as bp
from .base import BiopythonOverlayClass


class Atom(bp.Atom.Atom, BiopythonOverlayClass):
    child_class = None

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

class AtomDisorderException(Exception):
    pass


class BIAtom(BiopythonOverlayClass):
    child_class = None
    def __init__(self, data):
        if len(data) == 1:
            data = data[0]
        for k, v in data.items():
            if v == ".":
                data[k] = None

        essential_labs = [
            "group_PDB",
            "id",
            "type_symbol",
            "label_alt_id",
            "label_atom_id",
            "label_comp_id",
            "label_seq_id",
            "auth_asym_id",
            "auth_seq_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
        ]


        self.unused = {}
        for k, v in data.items():
            if not (k  in essential_labs):
                self.unused[k] = v
        self.misc = {}
        self._pdbx_PDB_ins_code = data["pdbx_PDB_ins_code"]
        self._read_id = int(data["id"])
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

        if "label_alt_id" not in data:
            data["label_alt_id"] = "."
        self.alt_id = data["label_alt_id"]
        if self.alt_id is None or self.alt_id == "None":
            self.alt_id = "." 

        if self.alt_id == ".":
            self.disordered = False
            self.doppelgangers = None
            self.favourite = None
        else:
            self.disordered = True
            self.doppelgangers = []
            self.favourite = True

    def __repr__(self):
        if self.disordered:
            if self.favourite:
                r = "\n<bi.{} id={}.{} b={} occupancy={} (disordred)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
                for datm in self.doppelgangers:
                    r += "\n - <bi.{} id={}.{} b={} occupancy={} (disordered)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
                return r
            else:
                return "<bi.{} id={}.{} b={} occupancy={} (disordred/not-favourite)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
        else:
            return "<bi.{} id={}> b={}".format(self.__class__.__name__, self.id, self.b)

    def __str__(self):
        return repr(self)

    def __iter__(self):
        if not self.disordered:
            self.i = None
            raise AtomDisorderException("Trying to iterate over a non-disordered atom")
        else:
            self.i = 0
            return self

    def __next__(self):
        if self.i > len(self.doppelgangers):
            self.i = None
            raise StopIteration
        else:
            if self.i == 0:
                self.i += 1
                return self
            else:
                self.i += 1
                return self.doppelgangers[self.i - 1]

    def set_bfactor(self, bfactor, doppelgangers=True):
        if self.disordered and self.doppelgangers is not None:
            targets = [a for a in self]
        else:
            targets = [self]

        for t in targets:
            t.b = float(bfactor)
        return self.b



    def set_misc(self, label, value):
        if self.disordered and self.doppelgangers is not None:
            targets = [a for a in self]
        else:
            targets = [self]

        for t in targets:
            if label in self.unused:
                t.unused[label] = value
            else:
                t.misc[label] = value
        return self.misc[label]


    def _mmcif_dict(self, include_unused=True, include_misc=False):
        if self.alt_id is None:
            self.alt_id = "."

        data = {
            "group_PDB": f"{self.type:6s}",
            "id": f"{self._read_id:4d}",
            "label_alt_id": f"{self.alt_id}",
            "label_seq_id": f"{self.atomnum:4d}",
            "type_symbol": f"{self.element:3s}",
            "label_atom_id": f"{self.name:3s}",
            "label_comp_id": f"{self.resname:4s}",
            "auth_seq_id": f"{self.resnum:4d}",
            "auth_asym_id": f"{self.chain:2s}",
            "Cartn_x": f"{self.x:6.3f}",
            "Cartn_y": f"{self.y:6.3f}",
            "Cartn_z": f"{self.z:6.3f}",
            "occupancy": f"{self.occupancy:6.3f}",
            "B_iso_or_equiv": f"{self.b:4.2f}",

        }
        if include_unused or include_misc:
            for k, v in self.unused.items():
                assert k not in data.keys()
                data[k] = f"{v}"

        if include_misc:
            for k, v in self.misc.items():
                data[k] = f"{v}"


        return data







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