import os, json
from ..utilities import *
from ..utilities import clean_string
from ..utilities.exceptions import *
from .mmcif import *
import numpy as np



class BIEntity(object):
    child_class = None
    extension = "structure"
    level = "structure"
    tmp_folder = "/tmp"

    def __init__(self, export_folder="./bioiain/exports", parent=None, use_tmp=False, **kwargs):
        from . import BIAtom, BIResidue, BIChain

        self.children = []
        self.paths = {
            "self": None, # This entity cif path
            "minimal": None, # This but only nice atoms (no headers or data)
            "parent": None, # Parent entity cif path
            "export_folder": export_folder, # Folder with all exports (default: "bioiain/exports")
            "top_folder": None, # Highest related folder
            "sub_folder": "", # Path of self under top_folder
        }
        self.data = {
            "info": {
                "code": None, # The code of this structure, if any
                "name": None, # The name of this structure (used mainly for file naming)
                "class": self.__class__.__name__,
            },
            "sequences": {
                "aa": None,
            },
            "symmetry": {
            }
        }
        self.headers = {
            "entry":{
                "id":None
            }
        }
        self.flags = {
            "loaded": False,
        }
        self.exporting = ["data", "paths", "flags"]

        #Properties
        self._com = None

        #CVectors
        self._cvectors = None

        # Children
        self._chains = None
        self._residues = None
        self._atoms = None
        self._mates = None

        # Crystal
        self._card = None
        self._parameters = None
        self._operations = None

        if parent is not None:
            self.paths["export_folder"] = parent.paths["export_folder"]
            self.paths["parent"] = parent.paths["self"]
            self.data["info"]["parent"] = repr(parent)
            self.headers = parent.headers

        if use_tmp:
            self.paths["export_folder"] = os.path.join(self.tmp_folder, self.paths["export_folder"] )




    def __repr__(self):
        if self.is_symmetry():
            return "<{}:{} id={} op={}>".format(self.__class__.__name__, self.code(), self.id(), self.op())
        return "<{}:{} id={}>".format(self.__class__.__name__, self.code(), self.id())

    def __str__(self):
        return repr(self)

    def __len__(self):
        return len(self.residues())

    def name(self):
        return self.data["info"]["name"]

    def set_name(self, name, append=False):
        if append and self.name() is not None:
            self.data["info"]["name"] = self.name() + f"_{name}"
        else:
            self.data["info"]["name"] = name

    def path(self, minimal=False):
        if not minimal:
            if self.paths.get("self", None) is None:
                self.export()
            return self.paths["self"]
        else:
            if self.paths.get("minimal", None) is None:
                self.export(minimal=True)
            return self.paths["minimal"]

    def folder(self):

        folders = []
        if self.paths["export_folder"] is not None:
            folders.append(self.paths["export_folder"])
        if self.paths["top_folder"] is not None:
            folders.append(self.paths["top_folder"])
        if self.paths["sub_folder"] is not None:
            folders.append(self.paths["sub_folder"])

        return os.path.join(*folders)

    def code(self):
        return self.data["info"]["code"]

    def id(self):
        return self.data["info"]["code"]

    def get_sequence(self, name="aa"):
        return self.data["sequences"][name]

    def set_sequence(self, name, seq):
        self.data["sequences"][name] = seq
        return self.data["sequences"][name]

    def has_flag(self, flag, value=None):
        if value is None:
            return flag in self.flags

        return self.flags.get(flag) == value

    def set_flag(self, flag, value=None):
        self.flags[flag] = value

    def sequence(self, force=False):
        if self.get_sequence() is None or force:
            seq = "".join([r.rn1 for r in self.residues()])
            self.data["sequences"]["aa"] = seq
        return self.get_sequence()

    def structure(self, code=None):
        from .structure import BIStructure
        if code is None:
            code = self.data["info"]["code"]
        return BIStructure.from_atoms(self._atoms, code, parent=self)

    def chains(self):
        return self.atoms(as_chains=True, hetatm=True)

    def residues(self):
        return self.atoms(ca_only=False, residues=True)

    def waters(self):
        return self.atoms(ca_only=False, water=True)

    def dna(self):
        return self.atoms(ca_only=False, dna=True)

    def ligands(self, relevant_only=True):
        return [l for l in self.atoms(ca_only=False, ligands=True) if l.relevant]

    def cvectors(self):
        if self._cvectors is None:
            self._calculate_cvectors()
        return self._cvectors

    def _calculate_cvectors(self):
        log(1, "Calculating CVectors for:", self.name())
        from ..aleph.vectors import CVector
        residues = self.residues()
        n_res = len(residues)
        cvector_list = []
        for n, res in enumerate(residues):
            if n == 0 or n == n_res -1:
                continue
            print(f"{n:4d}/{len(residues)-2:4d}", end="\r")
            cvector = CVector(residues[n-1], res, residues[n+1], params=self.params(), symops=self.symops(), entity_centre=self.com())
            cvector_list.append(cvector)
        log(2, f"n CVectors: {len(cvector_list)}")
        self._cvectors = cvector_list
        return cvector_list


    def atoms(self, ca_only=False, hetatm=False, ligands=False, residues=False, dna=False, water=False, force=False, group_by_residue=False, disordered=False, as_residues=False, chain=None, group_by_chain=False, as_chains=False, **kwargs):
        from .atom import _fix_disordered
        from .residue import build_res

        target_entities = []
        if as_residues:
            residues = True

        if residues:
            target_entities.append("residue")

        if dna:
            target_entities.append("dna")


        if water or ligands:
            hetatm = True
            ca_only = False

        if water:
            target_entities.append("water")

        if ligands:
            target_entities.append("ligand")


        atoms = self.all_atoms()

        if not disordered:
            atoms = _fix_disordered(atoms)

        if not hetatm:
            atoms = [a for a in atoms if a.type == "ATOM"]
        if ca_only:
            atoms = [a for a in atoms if a.name == "CA"]
        if chain is not None:
            atoms = [a for a in atoms if a.chain == chain]

        if as_chains or group_by_chain:
            from .chain import BIChain
            chain_list = {}
            for atom in atoms:
                if atom.complex in chain_list.keys():
                    chain_list[atom.complex].append(atom)
                else:
                    chain_list[atom.complex] = [atom]
            if as_chains:
                for ch, atms in chain_list.items():
                    chain_list[ch] = BIChain().from_atoms(atms, self.code(), ch, parent=self)
                chain_list = [ch for ch in chain_list.values()]
            return chain_list

        if group_by_residue or len(target_entities) > 0:
            atoms_by_res = {}
            for atom in atoms:
                if atom.id2[1:] in atoms_by_res:
                    atoms_by_res[atom.id2[1:]].append(atom)
                else:
                    atoms_by_res[atom.id2[1:]] = [atom]

            if len(target_entities) > 0:
                entities = []
                for resatms in atoms_by_res.values():
                    e = build_res(resatms)
                    if e.type in target_entities:
                        entities.append(e)
                return entities

            return atoms_by_res

        return atoms



    def remove_atom(self, atom):
        try:
            print(len(self._atoms))
            self._atoms.remove(atom)
            print(len(self._atoms))
            log(f"warning", f"Atom {atom} removed")
        except:
            log("warning", f"Atom {atom} not found, not removed")
            return False
        return True

    def set_symmetry(self):
        pass



    @classmethod
    def from_atoms(cls, atoms, code=None, share=True, **kwargs):
        self = cls(**kwargs)
        if share:
            self._atoms = atoms
        else:
            self._atoms = [a.copy() for a in atoms]
        if code is None and kwargs.get("parent", None) is not None:
            code = kwargs["parent"].code()
        if code is not None:
            self.data["info"]["code"] = clean_string(code).upper()
            self.data["info"]["code"] = code
        self.set_name(self.code())
        self.paths["top_folder"] = self.code()

        return self



    @classmethod
    def from_file(cls, filepath, code="auto", file_format="auto", force=False, **kwargs):
        log(1, "Loading from file:", filepath)
        self = cls(**kwargs)
        if self.has_flag("loaded", True):
            if force:
                log("warning: already loaded")
            else:
                raise exceptions.AlreadyLoaded()

        if file_format == "auto":
            file_format = filepath.split(".")[-1]

        if file_format not in ["cif", "pdb"]:
            raise exceptions.UnknownFormat(file_format)

        self._all_atoms(filepath=filepath, force=True, is_pdb=file_format == "pdb")

        if code == "auto":
            code = read_mmcif(filepath, subset="_entry")["_entry.id"]

        if code == "file" or code is None:
            code = filepath.split(".")[0]

        self.data["info"]["code"] = clean_string(code).upper()
        self.headers["entry"]["id"] = self.code()
        self.paths["top_folder"] = self.code()
        self.set_name(self.code())
        self.set_flag("fractional", False)
        self.set_flag("loaded", True)
        return self




    def from_biopython(self, entity):
        pass


    def com(self, force=False):
        from ..utilities.maths import find_com
        if self._com is None or force:
            self._com = find_com(self)
        return self._com


    def all_atoms(self):
        if self._atoms is None:
            self._all_atoms()
        return self._atoms

    def _all_atoms(self, filepath=None, force=False, is_pdb=False):
        from .atom import BIAtom


        if filepath is None:
            filepath = self.paths.get("self", None)

        if not hasattr(self, "_atoms"):
            force = True
        elif self._atoms is None:
            force = True

        if force:
            if filepath is None:
                self.export()
            if not os.path.exists(filepath):
                self.export()

            log(1, "Reading atoms from CIF:", filepath)
            mmcif= read_mmcif(filepath, subset=["_atom_site", "_cell", "_symmetry"])
            atoms=mmcif("_atom_site")
            self.headers["cell"] = mmcif("_cell")
            self.headers["symmetry"] = mmcif("_symmetry")
            self.data["symmetry"]["in_asu"] = True
            self._calculate_crystal()
            atoms = [BIAtom(a) for a in atoms]

            self._atoms = atoms
        return self._atoms

    def fix_headers(self):
        if self.headers["symmetry"].get("space_group_name_H-M", None) is not None:
            self.headers["symmetry"]["space_group_name_H-M"] = f"\'{self.headers["symmetry"]["space_group_name_H-M"]}\'"


    def export(self, minimal=False, cleanup=False, as_pdb=False, target_folder=None):

        if target_folder is None:
            target_folder = self.paths["export_folder"]
            custom_folder = False
        else:
            custom_folder = True

        fname = f"{self.name()}.{self.extension}"
        if minimal:
            fname += ".minimal"
        try:
            base_folder = os.path.join(target_folder, self.paths.get("top_folder", self.code()), self.paths["sub_folder"])
        except TypeError:
            print(self.paths)
            raise
        os.makedirs(base_folder, exist_ok=True)
        base_path = os.path.join(base_folder, fname)
        log(2, f"Exporting: {self} to {base_path}")
        if self.has_flag("is_fractional", True):
            log("Warning", "A fractional entity was about to be exported!")
            log("Warning", "An orthogonal copy was made for you and exported instead! (only for atoms)")
            orth = self.copy()._to_orthogonal()
        else:
            orth = self

        if minimal:
            minimal_path= orth._export_structure(base_path, headers=False, misc_fields=True, cleanup=True, as_pdb=as_pdb)
            if not custom_folder:
                self.paths["minimal"] = minimal_path
            return minimal_path
        else:
            path = orth._export_structure(base_path, headers=True, misc_fields=True, cleanup=cleanup, as_pdb=as_pdb)
            if not custom_folder:
                self.paths["self"] = path
                if not as_pdb:
                    self.set_flag("exported", True)
                    self.paths["data"] = self._export_data(base_path)
            return path

    def _export_data(self, filepath, mode="w") -> str:
        if filepath.endswith(".cif"):
            filepath = filepath.replace(".cif", ".json")
        elif not filepath.endswith(".json"):
            filepath += ".json"
        exp = {e: self.__getattribute__(e) for e in self.exporting if hasattr(self, e)}
        with open(filepath, mode) as f:
            f.write(json.dumps(exp, indent=4))
        return filepath


    def _export_structure(self, filepath:str, atoms:list=None, headers:bool=None, misc_fields:bool=True, cleanup=False, as_pdb=False) -> str:
        mode = "w"
        if atoms is None:
            if cleanup:
                atoms = self.atoms()
            else:
                atoms = self._all_atoms()
        if as_pdb:
            return write_pdb_atoms(atoms, filepath, mode=mode, end=True)
        else:
            if headers:
                for e in self.exporting:
                    d = getattr(self, e)
                    if e == "data":
                        for k, v in d.items():
                            write_dict(v, file_path=filepath, label=f"bi_{e}_{k}", mode=mode, name=self.name())
                            mode = "a"
                    else:
                        write_dict(d, file_path=filepath, label=f"bi_{e}", mode=mode, name=self.name())
                        mode = "a"
                for k, d in self.headers.items():
                    write_dict(d, file_path=filepath, label=k, mode=mode, name=self.name())
                    mode = "a"

            return write_atoms(atoms, filepath, name=self.name(), include_misc=misc_fields, mode=mode)

    @classmethod
    def recover_from_id(cls, code, endswith=None, full_name=None, **kwargs):
        placeholder = cls(**kwargs)
        path = os.path.join(placeholder.paths["export_folder"], code, placeholder.paths["sub_folder"])
        if full_name is not None:
            path = os.path.join(path, full_name)
        else:
            if not os.path.exists(path):
                raise StructureNotFound(path)
            for file in os.listdir(path):
                ext = file.split(".")[-1] 
                if ext != "json":
                    continue
                extension= file.split(".")[-2]
                if extension != placeholder.extension:
                    continue
                if endswith is not None:
                    if not file.split(".")[0].endswith(endswith):
                        continue
                path = os.path.join(path, file)
                return cls.recover_from_path(path, **kwargs)

        raise StructureNotFound(path)




    @classmethod
    def recover_from_path(cls, path, **kwargs):

        if path.endswith(".json"):
            data_path = path
            cif_path = data_path.replace(".json", ".cif")
        elif path.endswith(".cif"):
            cif_path = path
            data_path = cif_path.replace(".cif", ".json")
        else:
            cif_path = path+".cif"
            data_path = path+".json"

        raw = None
        if os.path.exists(cif_path):
            self = cls.from_file(cif_path, **kwargs)
        elif os.path.exists(data_path):
            raw = json.load(open(path, "r"))
            self = cls.from_file(raw["paths"]["self"], **kwargs)
        else:
            raise FileNotFoundError(path, cif_path, data_path)

        if os.path.exists(data_path) and raw is None:
            raw = json.load(open(path, "r"))
            for k, v in raw.items():
                setattr(self, k, v)
        else:
            log("Warning", "Data file not found for:", path)

        return self

    def _calculate_crystal(self):
        self._get_crystal_card()
        self._get_operations()

    def _get_crystal_card(self):
        cell = self.headers["cell"]

        a = float(cell["length_a"])
        b = float(cell["length_b"])
        c = float(cell["length_c"])
        alpha = float(cell["angle_alpha"])
        beta = float(cell["angle_beta"])
        gamma = float(cell["angle_gamma"])
        Z = float(cell["Z_PDB"])

        card = dict(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma, Z=Z)
        self._card = card
        return self._card

    def _get_operations(self):
        from ..utilities.space_groups import dictio_space_groups

        space_group_key = self.headers["symmetry"].get("Int_Tables_number")

        space_group_key = int(space_group_key)
        self._operations = dictio_space_groups[space_group_key]
        return self._operations

    def operations(self):
        if self._operations is None:
            self._get_operations()
        return self._operations

    def symops(self, n=None):
        if self._operations is None:
            self._get_operations()
        if n is None:
                return self._operations["symops"]
        else:
            return self._operations["symops"][n]

    def params(self):
        if self._parameters is None:
            self._calculate_parameters()
        return self._parameters


    def _calculate_parameters(self) -> dict:
        if self._card is None:
            self._get_crystal_card()

        card = self._card

        parameters = {}
        parameters["A"] = A = float(card["a"])
        parameters["B"] = B = float(card["b"])
        parameters["C"] = C = float(card["c"])
        parameters["alphaDeg"] = alphaDeg = float(card["alpha"])
        parameters["betaDeg"] = betaDeg = float(card["beta"])
        parameters["gammaDeg"] = gammaDeg = float(card["gamma"])
        parameters["alpha"] = alpha = (alphaDeg * 2 * np.pi) / 360
        parameters["beta"] = beta = (betaDeg * 2 * np.pi) / 360
        parameters["gamma"] = gamma = (gammaDeg * 2 * np.pi) / 360
        parameters["c_a"] = c_a = np.cos(alpha)
        parameters["c_b"] = c_b = np.cos(beta)
        parameters["c_g"] = c_g = np.cos(gamma)
        parameters["s_g"] = s_g = np.sin(gamma)
        parameters["q"] = q = np.sqrt(1 + 2 * c_a * c_b * c_g - c_a ** 2 - c_b ** 2 - c_g ** 2)
        parameters["uu"] = uu = s_g / (q * C)
        parameters["vv"] = vv = (c_b * c_g - c_a) / (q * B * s_g)
        parameters["uuy"] = uuy = 1 / (B * s_g)
        parameters["vvz"] = vvz = -1 * (c_g / (A * s_g))
        parameters["uuz"] = uuz = (c_a * c_g - c_b) / (q * A * s_g)
        parameters["vvy"] = vvy = 1 / A

        self._parameters = parameters
        return self._parameters

    def fragment(self, in_place=False, force=False):
        from ..aleph.fragments import FragmentedStructure
        if isinstance(self, FragmentedStructure):
            if not in_place:
                frag = self.copy()
            else:
                frag = self
        else:
            frag = FragmentedStructure.from_atoms(self.all_atoms(), parent=self, share=in_place, export_folder=self.paths["export_folder"])


        frag.fragments(force=force)
        return frag


    def copy(self):
        from copy import deepcopy
        new = self.__class__.from_atoms(self.all_atoms(), parent=self, share=False)
        new.data = deepcopy(self.data)
        new.flags = deepcopy(self.flags)
        new.set_flag("is_copy", True)
        new.set_flag("exported", False)
        return new

    def displace(self, distance:float|int|list[float|int]|tuple[float|int], inplace=True):
        if not inplace:
            self = self.copy()
        for a in self.all_atoms():
            a + distance
        return self

    def _to_fractional(self):
        if self.has_flag("is_fractional", True):
            raise AlreadyFractional(self)
        for a in self.all_atoms():
            a.to_frac(self.params())
        self.set_flag("is_fractional", True)

    def _to_orthogonal(self):
        if self.has_flag("is_fractional", False):
            raise AlreadyOrthogonal(self)
        for a in self.all_atoms():
            a.to_orth(self.params())
        self.set_flag("is_fractional", False)



    def _symmetry_operation(self, symop):
        was_orth = False
        if self.has_flag("is_fractional", False):
            was_orth = True
            self._to_fractional()

        for a in self.all_atoms():
            a.symop(self.symops(symop), self.params())

        if was_orth:
            self._to_orthogonal()
        return self

    def is_symmetry(self):
        return self.has_flag("is_symmetry", True)

    def op(self):
        return self.data["symmetry"].get("symop", None)


    def symmetry(self, symop, in_place=False):
        if not in_place:
            self = self.copy()

        self._symmetry_operation(symop)
        self.set_flag("is_symmetry", True)
        self.data["symmetry"]["in_asu"] = False
        self.data["symmetry"]["symop"] = symop
        return self


    def mates(self):
        mates = []
        for symop in self.symops().keys():
            mates.append(self.symmetry(symop))
        self._mates = mates
        return self._mates

    def show(self, execute=True, script=None):
        if script is None:
            from ..visualisation.pymol import PymolScript
            script = PymolScript(self.name())
        script.load(self.export(), self.name())
        script.spectrum(self.name())
        script.orient()
        script.write_script()
        if execute:
            script.execute()
        return script























