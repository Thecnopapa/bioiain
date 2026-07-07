from ..utilities.exceptions import *
from .mmcif import *
import numpy as np
from ..utilities import *


class BIEntity(object):
    child_class = None
    extension = "structure"
    level = "structure"
    tmp_folder = "/tmp"
    excluded_from_headers = ["_bi_*", "_atom_site", "_aleph_*","_cell", "_symmetry","_entry"]

    def __init__(self, export_folder=None, parent=None, use_tmp=False, model:str|int="*", **kwargs):
        if export_folder is None:
            export_folder = os.path.join(SUBDIR_NAME, "exports").strip()
        self.children = []
        self.paths = {
            "self": None, # This entity cif path
            "source": None,
            "minimal": None, # This but only nice atoms (no headers or data)
            "parent": None, # Parent entity cif path
            "export_folder": export_folder.strip(), # Folder with all exports (default: "bioiain/exports")
            "top_folder": None, # Highest related folder
            "sub_folder": "", # Path of self under top_folder
            "fasta": None, # Path to fasta
            "db": None, #mmseqs DB path
        }
        self.data = {
            "info": {
                "code": None, # The code of this structure, if any
                "name": None, # The name of this structure (used mainly for file naming)
                "class": self.__class__.__name__,
                "model": str(model)
            },
            "sequences": {
                "aa": None,
            },
            "symmetry": {},

            # Machine
            "embeddings": {},

            # Tools
            "SASA": {},
            "PISA":{},
            "DSSP":{},
        }
        self.headers = {
            "entry":{
                "id":None
            }
        }
        self.flags = {
            "loaded": False,
            "no_atoms": True,
        }
        self.exporting = ["data", "paths", "flags"]

        #Properties
        self._com = None
        self._kdtrees = {}


        #CVectors
        self._cvectors = None
        self._missing_cvectors = None

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
            self.paths["export_folder"] = parent.paths["export_folder"].strip()
            self.paths["parent"] = parent.paths["self"]
            self.data["info"]["parent"] = repr(parent)
            self.headers = parent.headers

        if use_tmp:
            self.paths["export_folder"] = os.path.join(self.tmp_folder, self.paths["export_folder"] )


    def clear_cahces(self):
        self._com = None
        self._kdtrees = {}
        self._chains = None
        self._residues = None
        self._atoms = None
        self._mates = None
        self._cvectors = None
        self._missing_cvectors = None

    def __repr__(self):
        if self.is_symmetry():
            return "<bi.{}:{} id={} op={} model={}>".format(self.__class__.__name__, self.code(), self.id(), self.op(), self.model())
        #return "<{}:{} id={} (len:{})>".format(self.__class__.__name__, self.code(), self.id(), len(self))
        return "<bi.{}:{} id={} model={}>".format(self.__class__.__name__, self.code(), self.id(), self.model())

    def __str__(self):
        return repr(self)

    def __len__(self):
        if self.has_flag("no_atoms", True):
            return 0
        return len(self.residues())

    def model(self):
        return self.data["info"]["model"]

    def avail_models(self):
        return list(set([str(a.model) for a in self._all_atoms()]))


    def set_model(self, model_num:str|int):
        model_num = str(model_num)
        assert model_num in self.avail_models(), f"Requested model({model_num}) does not exist in this entity: {self.avail_models()}"
        self.data["info"]["model"] = model_num
        return self

    def models(self, in_place=True):
        if in_place:
            for m in self.avail_models():
                self.set_model(m)
                yield self
        else:
            models = []
            for m in self.avail_models():
                if m == self.model():
                    models.append(self)
                else:
                    new_model = self.copy()
                    new_model.set_model(m)
                    models.append(new_model)
            return models

    def full_id(self):
        return f"{self.code()}"

    def name(self):
        return str(self.data["info"]["name"])

    def set_name(self, name, append=False):
        if append and self.name() is not None:
            self.data["info"]["name"] = self.name() + f"_{name}"
        else:
            self.data["info"]["name"] = str(name)

    def set_misc(self, key, value, force=True):
        for a in self.all_atoms():
            if not force:
                try:
                    a.get_misc(key)
                    continue
                except KeyError:
                    pass
            a.set_misc(key, value)

    def path(self, minimal=False, source=False):
        assert not (minimal and source), f"Trying to get both minimal and source paths at the same time"
        if minimal:
            if self.paths.get("minimal", None) is None:
                self.export(minimal=True)
            return self.paths["minimal"].strip()
        elif source:
            return self.paths["source"].strip()
        else:
            if self.paths.get("self", None) is None:
                self.export()
            return self.paths["self"].strip()

    def folder(self):

        folders = []
        if self.paths["export_folder"] is not None:
            folders.append(self.paths["export_folder"])
        if self.paths["top_folder"] is not None:
            folders.append(self.paths["top_folder"])
        if self.paths["sub_folder"] is not None:
            folders.append(self.paths["sub_folder"])

        return os.path.join(*folders).strip()

    def code(self):
        return str(self.data["info"]["code"])

    def id(self):
        return str(self.data["info"]["code"])

    def get_sequence(self, name=None):
        if name is None:
            return self.sequence()
        return self.data["sequences"].get(name, None)

    def set_sequence(self, name, seq):
        self.data["sequences"][name] = seq
        return self.data["sequences"][name]

    def has_flag(self, flag, value=None):
        if value is None:
            return flag in self.flags

        return self.flags.get(flag) == value

    def set_flag(self, flag, value=True):
        self.flags[flag] = value

    def sequence(self, force=False):
        if self.get_sequence("aa") is None or force:
            seq = "".join([r.rn1 for r in self.residues()])
            self.data["sequences"]["aa"] = seq
        return self.get_sequence("aa")

    def structure(self, code=None):
        from .structure import BIStructure
        if code is None:
            code = str(self.data["info"]["code"])
        return BIStructure.from_atoms(self._atoms, code, parent=self)

    def chains(self, sele:list|str=None, by_complex=False, **kwargs):
        return self.atoms(as_chains=True, hetatm=True, chain_sele=sele, by_complex=by_complex, **kwargs)

    def residues(self, need_main=True, need_backbone=True, **kwargs) -> list:
        return self.atoms(ca_only=False, residues=True, hetatm=True, residue_params={"need_main":need_main, "need_backbone": need_backbone}, **kwargs)

    def waters(self, **kwargs):
        return self.atoms(ca_only=False, water=True, **kwargs)

    def dna(self, **kwargs):
        return self.atoms(ca_only=False, dna=True, **kwargs)

    def ligands(self, relevant_only=True, **kwargs):
        return [l for l in self.atoms(ca_only=False, ligands=True, **kwargs) if l.relevant]

    def cvectors(self, vc_mode=None, return_missing=False, force=False):
        if self._cvectors is None or (vc_mode != getattr(self, "_vc_mode", None)) or force:
            self._calculate_cvectors(vc_mode=vc_mode)
        else:
            log(1, f"Using previously saved CVectors ({vc_mode})...")
        if return_missing:
            return self._cvectors, self._missing_cvectors
        return self._cvectors

    def _calculate_cvectors(self, vc_mode=None):
        log(1, "Calculating CVectors for:", self.name(), f"({vc_mode})")
        from ..aleph.vectors import CVector
        residues = self.residues(hetatm=True)
        n_res = len(residues)
        cvector_list = []
        missing_cvectors = []
        for n, res in enumerate(residues):
            if n == 0 or n == n_res -1:
                missing_cvectors.append(n)
                continue
            print(f"{n:4d}/{len(residues)-2:4d}", end="\r")
            cvector = CVector(residues[n-1], res, residues[n+1], params=self.params(), symops=self.symops(), entity_centre=self.com(), vc_mode=vc_mode)
            if cvector.trash:
                missing_cvectors.append(n)
                print()
                print(cvector)
                continue
            cvector_list.append(cvector)
        log(2, f"n CVectors: {len(cvector_list)}, trash={len(missing_cvectors)} ")
        self._missing_cvectors = missing_cvectors
        self._cvectors = cvector_list
        self._vc_mode = vc_mode
        return self._cvectors


    def atoms(self,
              ca_only=False,
              hetatm=False,
              ligands=False,
              residues=False,
              dna=False,
              water=False,
              hydrogens=False,
              group_by_residue=False,
              disordered=False,
              as_residues=False,
              chain=None,
              group_by_chain=False,
              as_chains=False,
              by_complex=False,
              model:int|str|None = None,
              verbose=False,
              match_all=False,
              residue_params={},
              **kwargs) -> list:
        from .atom import _fix_disordered
        from .residue import build_res

        target_entities = []
        if as_residues:
            residues = True

        if residues:
            target_entities.append("residue")

        if dna:
            target_entities.append("nucleotide")


        if water or ligands:
            hetatm = True
            ca_only = False

        if water:
            target_entities.append("water")

        if ligands:
            target_entities.append("ligand")


        atoms = self.all_atoms(verbose=verbose)

        if model is None:
            model = self.model()

        if not disordered:
            atoms = _fix_disordered(atoms)

        if (model is not None) and model not in ["*" "-"]:
            atoms = [a for a in atoms if str(a.model) == str(model)]

        if not hetatm:
            atoms = [a for a in atoms if a.type == "ATOM"]
        if not hydrogens:
            atoms = [a for a in atoms if a.element != "H"]
        if ca_only:
            atoms = [a for a in atoms if a.name == "CA"]
        if not water:
            atoms = [a for a in atoms if a.name != "HOH"]
        if chain is not None:
            atoms = [a for a in atoms if a.chain == chain]

        if as_chains or group_by_chain:
            chain_sele = kwargs.get("chain_sele", None)
            if chain_sele in ["*", "", "-"]:
                chain_sele = None
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
                #[print(ch.id(), type(ch.id())) for ch in chain_list.values()]
                #print((ch.id() if not by_complex else ch.complex()) for ch in chain_list.values())
                chain_list = list(chain_list.values())
                chain_list = [ch for ch in chain_list if len(ch.residues()) > 0]
                if chain_sele is not None:
                    if by_complex:
                        chain_list = [ch for ch in chain_list if ch.complex() in chain_sele]
                    else:
                        chain_list = [ch for ch in chain_list if ch.id() in chain_sele]

            return chain_list

        if group_by_residue or len(target_entities) > 0:
            atoms_by_res = {}
            for atom in atoms:
                if atom.id2()[1:-1] in atoms_by_res:
                    atoms_by_res[atom.id2()[1:-1]].append(atom)
                else:
                    atoms_by_res[atom.id2()[1:-1]] = [atom]

            if len(target_entities) > 0:
                entities = []
                for k, resatms in atoms_by_res.items():

                    try:
                        ent = build_res(resatms, parent=self)
                    except NoMatchingClass:
                        if verbose:
                            log("warning", "No matching class for:", k)
                            [log("warning", ra) for ra in resatms]
                        if match_all:
                            raise
                        else: 
                            continue
                    except NoMainAtomFound:
                        if verbose:
                            log("warning", "No main atom for:", k)
                            [log("warning", ra) for ra in resatms]
                        if residue_params.get("need_main", True):
                            continue

                    except NoBackbone as e:
                        if verbose:
                            log("warning", "No backbone for:", k)
                            [log("warning", ra) for ra in resatms]
                        if residue_params.get("need_backbone", True):
                            continue
                        ent = e.residue
                        
                    except Exception as e:
                        log("warning", e.__class__.__name__, e)
                        raise e

                    if getattr(ent, "type", None) in target_entities:
                        entities.append(ent)
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
            #self.data["info"]["code"] = str(code)
        self.set_name(self.code())
        self.paths["top_folder"] = self.code()
        self.set_flag("no_atoms", False)

        return self



    @classmethod
    def from_file(cls, filepath, code="auto", file_format="auto", force=False, check_existing=True, source=None, export=False, no_atoms=False, verbose=False, **kwargs):
        if verbose:
            if check_existing:
                log(1, "Loading from file:", filepath, f"(atoms={not no_atoms})")
            else:
                log(2, "Loading from file:", filepath, f"(atoms={not no_atoms})")
        if not os.path.exists(filepath):
            raise FileNotFoundError(filepath)
        self = cls(**kwargs)

        if file_format == "auto":
            file_format = filepath.split(".")[-1]

        if file_format not in ["cif", "pdb"]:
            raise UnknownFormat(file_format)
        if not no_atoms:
            self.set_flag("no_atoms", False)
            try:
                self._all_atoms(filepath=filepath, force=True, is_pdb=file_format == "pdb", **kwargs)
            except (StructureLoadException, CrystalError) as e:
                #log("Error", f"Structure not loaded: {filepath}", e.__class__.__name__)
                raise StructureLoadException(f"Structure not loaded: {filepath}", e.__class__.__name__)
                return None
        else:
            self.set_flag("no_atoms", True)

        if code == "auto":
            code = str(read_mmcif(filepath, subset="_entry")["_entry.id"])

        if code == "file" or code is None:
            code = filepath.split(".")[0]

        self.data["info"]["code"] = clean_string(code).upper()
        self.headers["entry"]["id"] = self.code()
        self.paths["top_folder"] = self.code()
        self.set_name(self.code())
        self.set_flag("fractional", False)

        if verbose:
            log(2, "Reading CIF data...")
        data = read_mmcif(filepath, subset="_bi_*", verbose=verbose)
        if len(data) > 0:
            #print(data)
            for k in data.keys():
                kk = k.split("_")[2]
                ss = "_".join(k.split("_")[3:])
                old_data = getattr(self, kk, {})
                new_data = old_data
                if ss is None or ss == "":
                    new_data = new_data | data[k]
                else:
                    new_data[ss] = new_data.get(ss, {}) | data[k]
                setattr(self, kk, new_data)
            self.set_flag("data_recovered", True)
        else:
            self.set_flag("data_recovered", False)

        if self.paths["self"] != filepath:
            self.paths["source"] = filepath
        if self.paths["sub_folder"] is None:
            self.paths["sub_folder"] = ""
        elif source is not None:
            self.paths["source"] = source


        if check_existing and not force:
            if verbose:
                log(2, "Checking previous exported data...")
            prev_path = self.export(dry=True)
            if os.path.exists(prev_path):
                if verbose:
                    log(3, "Recovering previously exported file...")
                try:
                    recovered_self = cls.from_file(prev_path, check_existing=False, source=filepath, no_atoms=no_atoms, **kwargs)
                    if recovered_self is None:
                        raise StructureRecoverException()
                    return recovered_self
                except Exception as e:
                    if verbose:
                        log("warning", e.__class__.__name__, e)
                        log("Error", "Recovery failed (returning new)")
                    raise e

        self._recover_cvectors()

        self.set_flag("loaded", True)
        self.set_flag("no_atoms", no_atoms)

        if export:
            self.export()
        return self


    def _recover_cvectors(self):
        pass


    def from_biopython(self, entity):
        pass


    def com(self, force=False):
        from ..utilities.maths import find_com
        if self._com is None or force:
            self._com = find_com(self)
        return self._com


    def all_atoms(self, *args, **kwargs):
        if self._atoms is None:
            self._all_atoms()
        return self._atoms

    def _all_atoms(self, filepath=None, force=False, require_crystal=True, verbose=False, **kwargs):
        from .atom import BIAtom


        if filepath is None:
            filepath = self.paths.get("self", None)

        if filepath is None:
            filepath = self.paths.get("source", None)
        #print(filepath)
        if not hasattr(self, "_atoms"):
            force = True
        elif self._atoms is None:
            force = True

        if force:
            if filepath is None:
                filepath = self.export()
            if not os.path.exists(filepath) and self.has_flag("no_atoms", False):
                filepath = self.export()
            filepath = filepath.strip()
            #print(filepath)
            if verbose:
                log(2, "Reading atoms from CIF:", filepath)
            try:
                mmcif= read_mmcif(filepath, subset=["_atom_site", "_cell", "_symmetry"])
            except:
                raise MissingCrystalInfo()
            atoms=mmcif("_atom_site")
            if atoms is None:
                raise StructureLoadException(f"No atoms in: {filepath}")
            self.headers["cell"] = mmcif("_cell")
            self.headers["symmetry"] = mmcif("_symmetry")
            self.data["symmetry"]["in_asu"] = True
            try:
                self._calculate_crystal()
            except CrystalError as e:
                if require_crystal:
                    raise e
            atoms = [BIAtom(a) for a in atoms]

            self._atoms = atoms
        self.set_flag("no_atoms", False)
        return self._atoms

    def fix_headers(self):
        if self.headers["symmetry"].get("space_group_name_H-M", None) is not None:
            self.headers["symmetry"]["space_group_name_H-M"] = f"\'{self.headers['symmetry']['space_group_name_H-M']}\'"


    def export(self, minimal=False, cleanup=False, as_pdb=False, target_folder=None, sufix=None, dry=False, all_headers=True, cvmatrix=True, cvectors=True, verbose=False):

        custom_folder = False
        if target_folder is None:
            target_folder = self.paths["export_folder"]
            target_folder = str(target_folder).strip()

        else:
            custom_folder = True
            target_folder = str(target_folder).strip()

        fname = str(self.name())
        if sufix is not None:
            custom_folder = True
            if sufix[0] in [".", "-", "_", "(",]:
                fname += sufix
            else:
                fname += "_"+sufix

        fname += f".{self.extension}"
        if minimal:
            fname += ".minimal"
        try:
            base_folder = os.path.join(target_folder, str(self.paths.get("top_folder", self.code())), self.paths["sub_folder"]).strip()
        except TypeError:
            print(self.paths)
            raise
        os.makedirs(base_folder, exist_ok=True)
        base_path = os.path.join(base_folder, fname)
        if dry:
            if verbose:
                log(3, "Calculating export path...")
            return base_path+".cif"
        if verbose:
            log(2, f"Exporting: {self} to {base_path}")
        if self.has_flag("is_fractional", True):
            log("Warning", "A fractional entity was about to be exported!")
            log("Warning", "An orthogonal copy was made for you and exported instead! (only for atoms)")
            orth = self.copy()._to_orthogonal()
        else:
            orth = self

        if minimal:
            minimal_path= orth._export_structure(base_path, headers=False, all_headers=False, misc_fields=True, cleanup=True, as_pdb=as_pdb, cvectors=False, cvmatrix=False, verbose=verbose)
            if not custom_folder:
                self.paths["minimal"] = minimal_path
            return minimal_path
        else:
            path = orth._export_structure(base_path, headers=True, all_headers=all_headers, misc_fields=True, cleanup=cleanup, as_pdb=as_pdb, cvmatrix=cvmatrix, cvectors=cvectors, verbose=verbose)
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


    def _export_structure(self, filepath:str, atoms:list=None, headers:bool=None, misc_fields:bool=True, cleanup=True, as_pdb=False, all_headers=True, cvectors=True, cvmatrix=True, verbose=False) -> str:
        if verbose:
            log(2,"Exporting structure...")
        mode = "w"
        if atoms is None:
            if cleanup:
                atoms = self.atoms()
            else:
                #log("Warning", "Exporting all atoms and misc fields might corrupt the file (cleanup=True recommended)")
                atoms = self._all_atoms()
        filepath = filepath.strip()
        if as_pdb:
            return write_pdb_atoms(atoms, filepath.strip(), mode=mode, end=True)
        else:

            full_headers = {}
            if all_headers and headers and self.paths.get("source", None) is not None:
                full_headers = read_mmcif(self.paths.get("source", None), exclude=self.excluded_from_headers).dict()
                #print(full_headers.keys())

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
                # print(type(self.headers))
                # print(self.headers)


                full_headers = full_headers | self.headers
                for k, d in full_headers.items():
                    if type(d) is list and len(d) > 1:
                        write_dict_list(d, file_path=filepath, label=k, mode=mode, name=self.name())
                        mode = "a"
                    else:
                        if type(d) is list:
                            try:
                                d = d[0]
                            except:
                                log("warning", "Trying to export empty list")
                                continue
                        write_dict(d, file_path=filepath, label=k, mode=mode, name=self.name())
                        mode = "a"

            #print(self)
            #log(3, "CVECTORS", cvectors, self._cvectors is not None)
            #print(self._cvectors)
            #log(3,"CVMATRIX", cvmatrix, getattr(self, "_cvmatrix", None) is not None)
            #print(getattr(self, "_cvmatrix", None))


            if cvectors and (self._cvectors is not None):
                if verbose:
                    log(3, "Exporting cvectors...")
                write_dict_list(self._cvectors, file_path=filepath, label="aleph_cvectors", mode=mode, name=self.name())
                mode = "a"
            if cvmatrix and (getattr(self, "_cvmatrix", None) is not None):
                if verbose:
                    log(3, "Exporting cvmatrix...")
                write_dict_list([v.closest_vp for v in self._cvmatrix.vectors], file_path=filepath, label="aleph_cvmatrix", mode=mode, name=self.name())
                mode = "a"
            if verbose:
                log(3, "Exporting atoms...")
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

    def write_fasta(self, seq_name=None, force=False):
        folder = os.path.join(self.paths["export_folder"], str(self.paths.get("top_folder", self.code())), self.paths["sub_folder"]).strip()
        if seq_name is None:
            fasta_path = os.path.join(folder, f"{self.name()}.fasta")
        else:
            fasta_path = os.path.join(folder, f"{self.name()}.{seq_name}.fasta")
        if os.path.exists(fasta_path) and not force:
            return fasta_path
        seq = self.get_sequence(name=seq_name)
        if seq is None:
            log("Warning", f"Sequence {seq_name} not found")
            raise SequenceNotFound()
        with open(fasta_path, "w") as f:
            f.write(f"> {self.name()}\n")
            f.write(seq)
            f.write("\n")

        if seq_name is None:
            self.paths["fasta"] = fasta_path
        return fasta_path

    def db(self, force=False):
        from ..utilities.sequences import MMSEQS2
        mmseqs = MMSEQS2(self.write_fasta(), folder=self.folder(), force=force)
        self.paths["db"] = mmseqs.db_path()
        return mmseqs



    # @classmethod
    # def recover_from_path(cls, cif_path, **kwargs):
    #
    #     # if path.endswith(".json"):
    #     #     data_path = path
    #     #     cif_path = data_path.replace(".json", ".cif")
    #     # elif path.endswith(".cif"):
    #     #     cif_path = path
    #     #     data_path = cif_path.replace(".cif", ".json")
    #     # else:
    #     #     cif_path = path+".cif"
    #     #     data_path = path+".json"
    #     #
    #     # if os.path.exists(cif_path):
    #     #     self = cls.from_file(cif_path, check_existing=False, **kwargs)
    #     #     if self.has_flag("data_recovered", True):
    #     #         log(3, "Recovering data from CIF...")
    #     # else:
    #     #     raise FileNotFoundError(path, cif_path)
    #     # elif os.path.exists(data_path):
    #     #     raw = json.load(open(path, "r"))
    #     #     self = cls.from_file(raw["paths"]["self"], check_existing=False, **kwargs)
    #     # else:
    #     #     raise FileNotFoundError(path, cif_path, data_path)
    #     #
    #     # if os.path.exists(data_path) and raw is None:
    #     #     log(3, "Recovering data from json...")
    #     #     raw = json.load(open(path, "r"))
    #     #     for k, v in raw.items():
    #     #         setattr(self, k, v)
    #     # elif raw is None:
    #     #     log("Warning", "Data file not found for:", path)
    #
    #     return self

    def _calculate_crystal(self):
        try:
            self._get_crystal_card()
            self._get_operations()
        except Exception as e:
            self.set_flag("crystal_error", True)
            raise CrystalError()

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

    def fragment(self, in_place=False, force=False, aleph_mode="original"):
        from ..aleph.fragments import FragmentedStructure
        if isinstance(self, FragmentedStructure):
            if not in_place:
                frag = self.copy()
            else:
                frag = self
        else:
            frag = FragmentedStructure.from_atoms(self.all_atoms(), parent=self, share=in_place, export_folder=self.paths["export_folder"], aleph_mode=aleph_mode)


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

    def av_sasa(self, force=False, **kwargs):
        if force or self.data["SASA"].get("average", None) is not None:
            av = 0
            total = 0
            sasas = self.sasa_list(**kwargs)
            for s in sasas:
                if s is not None:
                    av += s
                    total += 1
            if total == 0:
                av = None
            else:
                av /= total
            self.data["SASA"]["average"] = av
        return self.data["SASA"]["average"]

    def sasa_list(self, **kwargs):
        self.calculate_sasa(**kwargs)
        sasas = []
        for a in self.atoms(**kwargs):
            sasas.append(a.get_misc("SASA", None))
        return sasas

    def calculate_sasa(self, force=False, **kwargs):
        if not self.has_flag("sasa_calculated") or force:
            from ..tools.SASA import SASA
            sasa = SASA(**kwargs)
            sasa.compute(self, **kwargs)
        return self

    def assembly_level(self):
        pass

    def _calculate_pisa(self, **kwargs):
        from ..tools.PISA import PISA
        pisa = PISA(pisa_id=self.name(), **kwargs)


    def img3D(self, size=16, property:None|dict|list|str=None, mode="mean", distortion="none", plot=False, show_plot=False, embedding=None, residue_kwargs={}, as_embedding=False, gif=False):
        log(2, f"Generating 3D voxels...")
        log(3, f"Size: {size}x{size}x{size} ({size**3})")
        log(3, f"Distortion: {distortion}")
        assert size % 2 == 0

        if gif or show_plot:
            plot = True

        output = None
        residues = self.residues(**residue_kwargs)



        if type(property) is str:
            properties = [{"property":property}]
        elif type(property) is list:
            properties = [{"property":p} for p in property]
        elif property is None:
            properties = [{"property":None}]
        else:
            raise NotImplementedError()

        if embedding is not None:
            as_embedding = True
            import torch
            #print(embedding.size())
            embedding = embedding.resize(len(residues), embedding.size()[-1])
            log(3, embedding.size())
            assert len(residues) == embedding.size()[0], f"N residues({len(residues)}) does not match embedding size ({embedding.size()[0]})"

            properties.extend([{"property": f"emb_{n}", "is_tensor":True, "n": n} for n in range(embedding.size()[-1])])

        pixels = np.arange(size**3)
        #print(pixels.shape)
        count_grid = pixels.reshape(size, size, size)
        value_grid = pixels.reshape(size, size, size)
        #print(count_grid)

        count_grid.fill(0)
        value_grid.fill(0.0)
        value_grid = value_grid.astype(np.float32)


        #print(count_grid)
        coords = np.array([res.ca.coord for res in residues])
        x_coords = np.array([c[0] for c in coords])
        y_coords = np.array([c[1] for c in coords])
        z_coords = np.array([c[2] for c in coords])

        coord_sets = [x_coords, y_coords, z_coords]


        #print(coords.shape)
        x_range = x_coords.min(axis=0), x_coords.max(axis=0)
        y_range = y_coords.min(axis=0), y_coords.max(axis=0)
        z_range = z_coords.min(axis=0), z_coords.max(axis=0)

        ranges = [x_range, y_range, z_range]

        if distortion == "none":
            vertice_len = max(end - start for start, end in ranges)
        else:
            raise NotImplementedError(f"Distortion mode {distortion} not implemented")

        log(3, f"Cube size: {vertice_len:5.3f} cubic A")

        pixel_len = vertice_len / (size-2)
        log(3, f"Voxel size: {pixel_len:5.3f} cubic A")

        x_center, y_centre, z_centre = centres = [start + ((end-start)/2) for start, end in ranges]

        centered_coords = []
        for coord_set, centre in zip(coord_sets, centres):
            centered = coord_set - centre
            #print(centered.shape)
            centered_coords.append(centered)

        pixeled_coords = []
        for coord_set in centered_coords:
            pixeled = (coord_set // pixel_len ) + size/2
            #print(pixeled.shape)
            pixeled_coords.append(pixeled)
        #print(pixeled_coords)
        pixeled_coords = np.array([(x,y,z) for x, y, z in zip(*pixeled_coords)]).astype(np.int64)
        #print(pixeled_coords.shape)



        for prop in properties:
            count_grid.fill(0)
            value_grid.fill(0.0)
            value_grid = value_grid.astype(np.float32)

            pp = prop.get("property")
            if pp is None:
                continue
            mm = prop.get("mode", mode)

            log(3, f"Property: {pp} (mode: {mm})", end="\r")


            #print(value_grid[6][6])

            for nn, (coord, res) in enumerate(zip(pixeled_coords, residues)):
                #print(coord)
                count_grid[coord[0], coord[1], coord[2]] = count_grid[coord[0], coord[1], coord[2]] +1
                if pp != "count":
                    if prop.get("is_tensor", False):
                        emb_n = prop.get("n")
                        p = float(embedding[nn, emb_n].item())
                        #print(nn, emb_n, p)
                    else:
                        if pp.lower() in ["b", "b_factor", "bfactor"]:
                            p = res.bfactor()
                        else:
                            p = res.get_misc(pp, None)
                        if p is None:
                            p = 0.
                        p = float(p)

                    if mode == "max":
                        value_grid[coord[0], coord[1], coord[2]] = max(value_grid[coord[0], coord[1], coord[2]], p)
                    elif mode == "min":
                        value_grid[coord[0], coord[1], coord[2]] = min(value_grid[coord[0], coord[1], coord[2]], p)
                    elif mode in ["sum", "mean"]:
                        value_grid[coord[0], coord[1], coord[2]] = value_grid[coord[0], coord[1], coord[2]] + p
                    else:
                        raise NotImplementedError(f"Mode {mode} not implemented")


            if pp == "count":
                value_grid = count_grid
                value_grid = value_grid.astype(np.float32)

            if property is not None and mode == "mean":
                #print(count_grid[:6, :6])
                #print(value_grid[:6, :6])
                #print("meaning")
                value_grid = np.divide(value_grid, count_grid, where=count_grid > 0.5) # Floating point error messing tings without this 0.5
                #print(value_grid[:6, :6])

            if as_embedding:
                import torch
                t = torch.tensor(value_grid)
                t = t.reshape(-1, *t.shape)
                #log(3,"Out:", t.shape)
            else:
                t = value_grid

            if output is None or not as_embedding:
                output = t
            else:
                output = torch.cat((output, t))
            #print(output.shape)


            if plot:
                gif_path = None
                if gif:
                    gif_path = os.path.join(self.folder(), "gifs", f"{self.full_id()}_{pp}")
                from src.bioiain.visualisation import voxels3d
                voxels3d(value_grid, count_grid, show_plot=show_plot, title=pp, gif_path=gif_path)



        print()
        return output
