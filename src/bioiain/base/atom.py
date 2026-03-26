import os, json

from ..utilities.exceptions import *
from ..utilities import *
import numpy as np


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
            "label_asym_id",
            "label_entity_id",
            "auth_atom_id",
            "auth_comp_id",
            "auth_seq_id",
            "auth_asym_id",
            "Cartn_x",
            "Cartn_y",
            "Cartn_z",
            "occupancy",
            "B_iso_or_equiv",
            "pdbx_PDB_model_num",

        ]


        self.misc = {}
        for k, v in data.items():
            if not (k  in essential_labs):
                self.misc[k] = v


        #ATOM
        self.atomnum = int(data["id"])
        self.type = data["group_PDB"] # ATOM / HETATM
        self.element = data["type_symbol"].strip() #C
        self.name = data["label_atom_id"].strip() # CA (Auto)
        self.name2 = data["auth_atom_id"].strip() # CA (Given)
        self.prime = False
        if "'" in self.name:
            self.prime = True
            self.name = self.name.replace("'", "").strip()
            self.name2 = self.name2.replace("'", "").strip()

        #RES
        self.resname = data["label_comp_id"] # Auto
        self.resname2 = data["auth_comp_id"] # Given
        self.resseq = data["label_seq_id"] # Auto
        if self.resseq is not None: self.resseq = int(self.resseq)
        self.resnum = data["auth_seq_id"] # Given
        if self.resnum is not None: self.resnum = int(self.resnum)

        #CHAIN
        self.chain = data["label_asym_id"] # Auto
        self.complex = data["auth_asym_id"] # Given
        self.entity = int(data["label_entity_id"])
        self.model = int(data["pdbx_PDB_model_num"])

        #PROPERTIES
        self.id = (self.name,  self.resnum, self.complex) # Ambiguous (label)
        self.id2 = (self.name,  self.resseq, self.chain) # Not ambiguous (auth)
        self.id3 = (self.atomnum, self.type, self.element, self.name, self.resseq, self.chain, self.model) # Unique
        self.x = float(data["Cartn_x"])
        self.y = float(data["Cartn_y"])
        self.z = float(data["Cartn_z"])
        self.occupancy = float(data["occupancy"])
        self.b = float(data["B_iso_or_equiv"])
        self.coord = (self.x, self.y, self.z)

        #CRYSTAL
        self.is_fractional = False

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
                    r += "\n - <bi.{} id={}.{} b={} occupancy={} (disordered)>".format(self.__class__.__name__, datm.id, datm.alt_id, datm.b, datm.occupancy)
                return r
            else:
                return "<bi.{} id={}.{} b={} occupancy={} (disordred/not-favourite)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
        else:
            return "<bi.{} id={} b={}>".format(self.__class__.__name__, self.id, self.b)

    def __repr__(self):
        if self.disordered:
            if self.favourite:
                return "<bi.{} id={}.{} b={} occupancy={} (disordred)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
            else:
                return "<bi.{} id={}.{} b={} occupancy={} (disordred/not-favourite)>".format(self.__class__.__name__, self.id, self.alt_id, self.b, self.occupancy)
        else:
            return "<bi.{} id={} b={}>".format(self.__class__.__name__, self.id, self.b)

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

    def set_coord(self, coord):
        for t in self:
            t.coord = float(coord)
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

    def __add__(self, other):
        if type(other) in (int, float):
            if len(other) == 1:
                a = (other, other, other)
            elif len(other) == 3:
                a = other
            else:
                raise TypeError
            self.x += a[0]
            self.y += a[1]
            self.z += a[2]
            self.coord = self.x, self.y, self.z
            return self
        else:
            raise NotImplementedError()

    def __sub__(self, other):
        if type(other) in (int, float):
            if len(other) == 1:
                a = (other, other, other)
            elif len(other) == 3:
                a = other
            else:
                raise TypeError
            self.x -= a[0]
            self.y -= a[1]
            self.z -= a[2]
            self.coord = self.x, self.y, self.z
            return self
        else:
            raise NotImplementedError()

    def copy(self):
        from copy import deepcopy
        new = deepcopy(self)
        return new

    @staticmethod
    def _to_frac(coord, params):
        x, y, z = coord
        nx = (x * params["vvy"]) + (y * params["vvz"]) + (z * params["uuz"])
        ny = (y * params["uuy"]) + (z * params["vv"])
        nz = z * params["uu"]
        return nx, ny, nz

    @staticmethod
    def _to_orth(coord, params):
        t1, t2, t3 = coord
        tz = t3 / params["uu"]
        ty = (t2 - tz * params["vv"]) / params["uuy"]
        tx = (t1 - ty * params["vvz"] - tz * params["uuz"]) / params["vvy"]
        return tx, ty, tz

    def to_frac(self, params):
        if self.is_fractional:
            log("Warining", f"Atom: {self} is already fractional")
            raise AlreadyFractional(self)

        self._orth_coords = self.coord

        nx, ny, nz = self._to_frac(self.coord, params)

        self.x = nx
        self.y = ny
        self.z = nz
        self.coord = (self.x, self.y, self.z)
        self.is_fractional = True

        return self

    def to_orth(self, params):
        if not self.is_fractional:
            log("Warining", f"Atom: {self} is already orthogonal")
            raise AlreadyOrthogonal(self)


        self._frac_coords = self.coord

        tx, ty, tz = self._to_orth(self.coord, params)

        self.x = tx
        self.y = ty
        self.z = tz
        self.coord = (self.x, self.y, self.z)
        self.is_fractional = False

        return self

    def at(self, symop, params, centre=None, centre_is_frac=False):
        from copy import deepcopy

        if not self.is_fractional:
            p1 = np.array(self._to_frac(self.coord, params))
            p2 = np.array(self._to_frac(self.coord, params))
            was_orth = True
        else:
            was_orth = False
            p1 = np.array(self.coord)
            p2 = np.array(self.coord)

        if centre is None:
            if self.is_fractional:
                centre = np.array(self.coord)
            else:
                centre = np.array(self._to_frac(self.coord, params))
        else:
            if centre_is_frac:
                centre = np.array(centre)
            else:
                centre = np.array(self._to_frac(centre, params))
                


        #print(self.coord)
        print("P1", p1)
        print("P2", p2)

        p1s = np.array(self._symop(p1, symop, params))
        p2s = np.array(self._symop(p2, symop, params)) + 99.5
        print("P1s", p1s)
        print("P2s", p2s)



        print("centre", centre)
        delta = ( ( p2s - np.array(centre) ) % 1 ) - 0.5
        print("delta", delta)
        new = np.array(centre) + delta
        print("P1s after", p1s)
        print("P2s after", p2s)
        print("new", new)
        print(new - p2s)
        position = (new - p2s + 99.5)
        print("position", position)
        for p in position:
            assert p % 1 == 0
        position = tuple([int(p) for p in position])



        #same = sum(abs(p1 - p2)) <= 0.0000001
        #print("same", same)
        same = False

        if was_orth:
            p1s = self._to_orth(p1s, params)
            if not same:
                new = self._to_orth(new, params)

        #p1 = tuple([float(c) for c in p1])
        if same:
            return p1s

        #p2 = tuple([float(c) for c in p2])
        return p1s, new


    def symop(self, symop, params):
        if not self.is_fractional:
            self.to_frac(params)
            was_orth = True
        else:
            was_orth = False

        nx, ny, nz = self._symop(self.coord, params)

        self.x = nx
        self.y = ny
        self.z = nz
        self.coord = (self.x, self.y, self.z)
        if was_orth:
            self.to_orth(params)

        return self


    def _symop(self, coord, symop, params):

        rot = symop["rot"]
        tra = symop["tra"]

        x, y, z = coord

        nx = (rot[0][0] * x) + (rot[0][1] * y) + (rot[0][2] * z) + tra[0]
        ny = (rot[1][0] * x) + (rot[1][1] * y) + (rot[1][2] * z) + tra[1]
        nz = (rot[2][0] * x) + (rot[2][1] * y) + (rot[2][2] * z) + tra[2]

        return nx, ny, nz


    @staticmethod
    def _none_point(val):
        if val is None:
            return "."
        return str(val)

    def _mmcif_dict(self, include_misc=True):

        try:
            data = {}
            data["group_PDB"] = f"{self.type:6s}"
            data["id"] = f"{self.atomnum:4d}"


            data["label_alt_id"] = f"{self._none_point(self.alt_id):>1s}"
            data["type_symbol"] = f"{self._none_point(self.element):<2s}"




            if self.prime and not self.name.endswith("'"): data["label_atom_id"] = f"\"{self.name:>2s}'\""
            else: data["label_atom_id"] = f"{self._none_point(self.name):<5s}"

            data["label_comp_id"] =f"{self._none_point(self.resname):<3s}"
            data["label_seq_id"] = f"{self._none_point(self.resseq):>3s}"
            data["label_asym_id"] = f"{self._none_point(self.chain):1s}"
            data["label_entity_id"] = f"{self._none_point(self.entity):1s}"

            if self.prime and not self.name2.endswith("'"): data["auth_atom_id"] = f"\"{self._none_point(self.name2):>2s}'\""
            else: data["auth_atom_id"] = f"{self._none_point(self.name2):<5s}"

            data["auth_comp_id"] =f"{self._none_point(self.resname2):<3s}"
            data["auth_seq_id"] = f"{self._none_point(self.resnum):>3s}"
            data["auth_asym_id"] = f"{self._none_point(self.complex):1s}"

            data["pdbx_PDB_model_num"] = f"{self._none_point(self.model):>2s}"
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
                data[k] = f"{self._none_point(v):>3s}"


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
