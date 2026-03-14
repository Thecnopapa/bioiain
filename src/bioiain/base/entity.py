import os, json
from ..utilities import *
from ..utilities import clean_string, d3to1
from .mmcif import *



class BIEntity(object):
    child_class = None
    extension = "structure"
    level = "structure"
    tmp_folder = "/tmp/bioiain"
    use_tmp = False

    def __init__(self, export_folder="bioiain/exports", parent=None, **kwargs):
        from . import BIAtom, BIResidue, BIChain

        self.children = []
        self.paths = {
            "self": None, # This entity cif path
            "parent": None, # Parent entity cif path
            "export_folder": export_folder, # Folder with all exports (default: "bioiain/exports")
            "top_folder": None, # Highest related folder
            "sub_folder": "", # Path of self under top_folder
        }
        self.data = {
            "info": {
                "code": None, # The code of this structure, if any
                "name": None, # The name of this structure (used mainly for file naming)
            },
            "sequences": {
                "aa": None,
            }
        }
        self.flags = {
            "loaded": False,
        }
        self.exporting = ["data", "paths", "flags"]
        self._chains = None
        self._residues = None
        self._atoms = None

        if parent is not None:
            self.paths["parent"] = parent.paths["self"]
            self.data["info"]["parent"] = repr(parent)


    def __repr__(self):
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

    def code(self):
        return self.data["info"]["code"]

    def id(self):
        return self.data["info"]["code"]

    def get_sequence(self, name="aa"):
        return self.data["sequences"][name]

    def has_flag(self, flag, value=None):
        if value is None:
            return flag in self.flags

        return self.flags.get(flag) == value

    def set_flag(self, flag, value=None):
        self.flags[flag] = value

    def sequence(self, force=False):
        if self.get_sequence() is None or force:
            seq = "".join([d3to1(r.resname) for r in self.residues()])
            self.data["sequences"]["aa"] = seq
        return self.get_sequence()

    def structure(self, code=None):
        from .structure import BIStructure
        if code is None:
            code = self.data["info"]["code"]
        return BIStructure.from_atoms(self._atoms, code)

    def chains(self):
        return self.atoms(as_chains=True, hetatm=True)

    def residues(self, force=False, chain=None, **kwargs):
        return self.atoms(ca_only=False, as_residues=True, force=force, chain=chain, **kwargs)




    def atoms(self, ca_only=False, hetatm=False, force=False, group_by_residue=False, disordered=False, as_residues=False, chain=None, group_by_chain=False, as_chains=False, **kwargs):
        from .atom import _fix_disordered
        from .residue import BIResidue
        atoms = self._all_atoms(force=force)

        if not disordered:
            atoms = _fix_disordered(atoms)

        if not hetatm:
            atoms = [a for a in atoms if a.type != "HETATM"]
        if ca_only:
            atoms = [a for a in atoms if a.name == "CA"]
        if chain is not None:
            atoms = [a for a in atoms if a.chain == chain]

        if as_chains or group_by_chain:
            from .chain import BIChain
            chain_list = {}
            for atom in atoms:
                if atom.chain in chain_list.keys():
                    chain_list[atom.chain].append(atom)
                else:
                    chain_list[atom.chain] = [atom]
            if as_chains:
                for ch, atms in chain_list.items():
                    chain_list[ch] = BIChain().from_atoms(atms, self.code(), ch, parent=self)
                chain_list = [ch for ch in chain_list.values()]
            return chain_list

        if group_by_residue or as_residues:
            atoms_by_res = {}
            for atom in atoms:
                if atom.id[1:] in atoms_by_res:
                    atoms_by_res[atom.id[1:]].append(atom)
                else:
                    atoms_by_res[atom.id[1:]] = [atom]


            atoms = atoms_by_res
            if as_residues:

                atoms = [BIResidue(resatms) for resatms in atoms.values()]

        return atoms

    def export(self):
        root = "."
        if self.use_tmp:
            root = self.tmp_folder
        fname = f"{self.name()}.{self.extension}"
        try:
            base_folder = os.path.join(root, self.paths["export_folder"], self.paths["top_folder"], self.paths["sub_folder"])
        except TypeError:
            print(self.paths)
            raise
        os.makedirs(base_folder, exist_ok=True)
        base_path = os.path.join(base_folder, fname)
        self.paths["self"] = self._export_structure(base_path)
        self.paths["data"] = self._export_data(base_path)
        return base_path






    @classmethod
    def from_atoms(cls, atoms, code=None, **kwargs):
        self = cls(**kwargs)
        self._atoms = atoms
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
        self.paths["top_folder"] = self.code()
        self.set_name(self.code())
        return self




    def from_biopython(self, entity):
        pass




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
                self._export()
            if not os.path.exists(filepath):
                self._export()

            log(1, "Reading atoms from CIF:", filepath)
            atoms= read_mmcif(filepath, subset=["_atom_site"])("_atom_site")
            atoms = [BIAtom(a) for a in atoms]

            self._atoms = atoms
        return self._atoms


    def _export_data(self, filepath, mode="w") -> str:
        if filepath.endswith(".cif"):
            filepath = filepath.replace(".cif", ".json")
        elif not filepath.endswith(".json"):
            filepath += ".json"
        exp = {e: self.__getattribute__(e) for e in self.exporting if hasattr(self, e)}
        with open(filepath, mode) as f:
            f.write(json.dumps(exp, indent=4))
        return filepath


    def _export_structure(self, filepath, atoms=None, headers=None, misc_fields=True) -> str:
        mode = "w"
        if atoms is None:
            atoms = self.atoms()
        if headers:
            #write_headers(self, filepath, mode=mode)
            mode = "a"
            pass
        return write_atoms(atoms, filepath, name=self.name(), include_misc=misc_fields, mode=mode)




    def recover(self, data_path):
        pass





