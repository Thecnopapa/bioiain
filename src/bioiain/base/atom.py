import os, json

from ..utilities.exceptions import *


class BIAtom(object):
    child_class = None
    def __init__(self, data):
        if len(data) == 1:
            data = data[0]
        for k, v in data.items():
            if v == "." or v == "?" or v=="None":
                #print(k, "is empty for:", data["id"])
                #print(data)
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


        self.misc = {}
        for k, v in data.items():
            if not (k  in essential_labs):
                self.misc[k] = v


        #ATOM
        self.atomnum = int(data["id"])
        self.type = data["group_PDB"] # ATOM / HETATM
        self.element = data["type_symbol"] #C
        self.name = data["label_atom_id"] # CA
        #RES
        self.resname = data["label_comp_id"]
        self.resseq = data["label_seq_id"] # Auto
        if self.resseq is not None: self.resseq = int(self.resseq)
        self.resnum = data["auth_seq_id"] # Given
        if self.resnum is not None: self.resnum = int(self.resnum)
        #CHAIN
        self.chain = data["auth_asym_id"]

        #PROPERTIES
        self.id = (self.name,  self.resnum, self.chain)
        self.id2 = (self.name,  self.resseq, self.chain)
        self.x = float(data["Cartn_x"])
        self.y = float(data["Cartn_y"])
        self.z = float(data["Cartn_z"])
        self.occupancy = float(data["occupancy"])
        self.b = float(data["B_iso_or_equiv"])
        self.coord = (self.x, self.y, self.z)

        #DISORDER
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



    def print_full(self):
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

    def __repr__(self):
        if self.disordered:
            if self.favourite:
                return "<bi.{} id={}.{} b={} occupancy={} (disordred)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
            else:
                return "<bi.{} id={}.{} b={} occupancy={} (disordred/not-favourite)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
        else:
            return "<bi.{} id={}> b={}".format(self.__class__.__name__, self.id, self.b)

    def __str__(self):
        return repr(self)

    def __iter__(self):
        if not self.disordered:
            self.i = 0
            #raise AtomDisorderException("Trying to iterate over a non-disordered atom")
            return self
        else:
            self.i = 0
            return self

    def __next__(self):
        if not self.disordered:
            if self.i == 0:
                self.i += 1
                return self
            else:
                self.i = None
                raise StopIteration

        if self.i > len(self.doppelgangers):
            self.i = None
            raise StopIteration
        else:
            if self.i == 0:
                self.i += 1
                return self
            else:
                self.i += 1
                #print(len(self.doppelgangers), self.i, self.i-2)
                return self.doppelgangers[self.i - 2]

    def set_bfactor(self, bfactor):
        for t in self:
            t.b = float(bfactor)
        return self.b



    def set_misc(self, label, value):

        for t in self:
            t.misc[label] = value
        return self.misc[label]

    def get_misc(self, label, keyerror="undefined"):
        if label in self.misc:
            return self.misc[label]
        else:
            if keyerror == "undefined":
                raise KeyError
            else:
                return keyerror

    def pdb_string(self, new_id=None):
        record_name = f"{self.type:<6s}"[-6:]
        if new_id is None:
            atom_serial_number = f"{self.atomnum:5d}"[-5:]
        else:
            atom_serial_number = f"{new_id:5d}"[-5:]
        if len(self.element) >= 2 or len(self.name) >= 4:
            atom_name = f"{self.name:<4s}"[-4:]
        else:
            atom_name = f" {self.name:<3s}"[-4:]
        if self.alt_id is None or self.alt_id==".":
            alternate_location_indicator = " "
        else:
            alternate_location_indicator = f"{self.alt_id:1s}"[0]
        residue_name = f"{self.resname:>3s}"[:3]
        chain_identifier = f"{self.chain:1s}"[0]
        residue_sequence_number = f"{self.resnum:>4d}"[-4:]
        code_for_insertions_of_residues = " " # ??
        xcoord = f"{self.x:>8.3f}"[:8]
        ycoord = f"{self.y:>8.3f}"[:8]
        zcoord = f"{self.z:>8.3f}"[:8]
        occupancy = f"{self.occupancy:>6.2f}"[-6:]
        temperature_factor = f"{self.b:>6.2f}"[-6:]
        segment_identifier = f"    " # ??
        element_symbol = f"{self.element:>2s}"[:2]
        charge = "  " # ??



        string = ""
        string+=record_name # 1-6
        string+=atom_serial_number # 7-11
        string+=" " # 12
        string+=atom_name # 13-16
        string+=alternate_location_indicator # 17
        string+=residue_name # 18-20
        string+=" " # 21
        string+=chain_identifier # 22
        string+=residue_sequence_number # 23-26
        string+=code_for_insertions_of_residues # 26
        string+="   " # 28-30
        string+=xcoord # 31-38
        string+=ycoord # 39-46
        string+=zcoord # 47-54
        string+=occupancy # 55-60
        string+=temperature_factor # 61-66
        string+= "      " # 67-72
        string+=segment_identifier # 73-76
        string+=element_symbol # 77-78
        string+=charge # 79-80

        return string


    def _mmcif_dict(self, include_misc=True):

        try:
            data = {}
            data["group_PDB"] = f"{self.type:6s}"
            data["id"] = f"{self.atomnum:4d}"
            if self.alt_id is None: data["label_alt_id"] = "."
            else: data["label_alt_id"] = f"{self.alt_id}"
            if self.resseq is None: data["label_seq_id"] = "."
            else: data["label_seq_id"] = f"{self.resseq:4d}"
            data["type_symbol"] = f"{self.element:3s}"
            data["label_atom_id"] = f"{self.name:3s}"
            data["label_comp_id"] =f"{self.resname:4s}"
            data["auth_seq_id"] = f"{self.resnum:4d}"
            data["auth_asym_id"] = f"{self.chain:2s}"
            data["Cartn_x"] = f"{self.x:7.3f}"
            data["Cartn_y"] = f"{self.y:7.3f}"
            data["Cartn_z"] = f"{self.z:7.3f}"
            data["occupancy"] = f"{self.occupancy:6.3f}"
            data["B_iso_or_equiv"] = f"{self.b:6.2f}"


        except Exception as e:
            print(self)
            print(self.__dict__)
            raise e


        if include_misc:
            for k, v in self.misc.items():
                if v is None: v = "."
                data[k] = f"{v}"


        return data





def _fix_disordered(atoms):
    fixed_atoms = []
    for atom in atoms:
        if atom.disordered:
            for a in fixed_atoms:
                if a.id2 == atom.id2:
                    if a.resseq is None and not (a.resnum == atom.resnum):
                        continue
                    try:
                        assert a.disordered
                    except:
                        print(">>>")
                        print(a, a.resseq, a.resnum)
                        print("###")
                        print(atom, atom.resseq, atom.resnum)
                        print("<<<")
                        raise
                    a.doppelgangers.append(atom)
                    a.favourite = True
                    atom.favourite = False
                    atom.doppelgangers = None
                    break
            if atom.favourite:
                fixed_atoms.append(atom)
        else:
            fixed_atoms.append(atom)
    return fixed_atoms
