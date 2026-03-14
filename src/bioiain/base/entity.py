import os, json
from ..utilities import *
from ..utilities import clean_string
from .mmcif import *



class BIEntity(object):
    child_class = None
    extension = "structure"
    level = "structure"
    tmp_folder = "/tmp/bioiain"
    use_tmp = False

    def __init__(self, export_folder="bioiain/exports"):
        from . import BIAtom, BIResidue, BIChain

        self.children = []
        self.paths = {
            "self": None, # This entity cif path
            "parent": None, # Parent entity cif path
            "export_folder": export_folder, # Folder with all exports (default: "bioiain/exports")
            "top_folder": None, # Highest related folder
            "sub_folder": None, # Path of self under top_folder
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
        self.id = None
        self._chains = None
        self._residues = None
        self._atoms = None


    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def name(self):
        return self.data["info"]["name"]

    def code(self):
        return self.data["info"]["code"]

    def sequence(self, name="aa"):
        return self.data["sequences"][name]

    def has_flag(self, flag, value=None):
        return self.flags.get(flag, None) == value

    def set_flag(self, flag, value=None):
        self.flags[flag] = value





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

    def atoms(self, ca_only=False, hetatm=False, force=False, group_by_residue=False, disordered=False, as_residues=False, chain=None, **kwargs):
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
        if group_by_residue or as_residues:
            atoms_by_res = {}
            residues_with_ca = []
            for atom in atoms:
                if atom.resnum in atoms_by_res:
                    atoms_by_res[atom.resnum].append(atom)
                else:
                    atoms_by_res[atom.resnum] = [atom]

                if atom.name == "CA":
                    residues_with_ca.append(atom.resnum)
            for key in list(atoms_by_res.keys()):
                try:
                    residues_with_ca.remove(key)
                except:
                    atoms_by_res.pop(key)

            atoms = atoms_by_res
            if as_residues:
                atoms = [BIResidue(resatms) for resatms in atoms.values()]

        return atoms


    def from_file(self, filepath, code="auto", level="structure", file_format="auto", force=False, export_folder=None):

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
        print(self.code())

        return self



    def read_pdb(self, filepath):
        pass

    def _from_biopython(self, entity):
        pass

    def export(self):
        root = "."
        if self.use_tmp:
            root = self.tmp_folder
        fname = f"{self.name()}.{self.extension}"
        base_folder = os.path.join(root, self.paths["export_folder"], self.paths["top_folder"], self.paths["sub_folder"])
        os.makedirs(base_folder, exist_ok=True)
        base_path = os.path.join(base_folder, fname)
        self.paths["self"] = self._export_structure(base_path)
        self.paths["data"] = self._export_data(base_path)
        return base_path

    def _export_data(self, filepath, mode="w") -> str:
        if filepath.endswith(".cif"):
            filepath = filepath.replace(".cif", ".json")
        elif not filepath.endswith(".json"):
            filepath += ".json"
        exp = {e: self.__getattribute__(e) for e in self.exporting if hasattr(self, e)}
        with open(filepath, mode) as f:
            f.write(json.dumps(exp, indent=4))
        return filepath


    def _export_structure(self, filepath, atoms=None, headers=None, unused_fields=False, misc_fields=False) -> str:
        mode = "w"
        if atoms is None:
            atoms = self.atoms()
        if headers:
            #write_headers(self, filepath, mode=mode)
            mode = "a"
            pass
        return write_atoms(self, filepath, include_unused=unused_fields, include_misc=misc_fields, mode=mode)




    def recover(self, data_path):
        pass





