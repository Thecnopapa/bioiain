import builtins
import os
from copy import deepcopy

import Bio.PDB as bp
import json
from ..utilities.logging import log
from typing_extensions import Self


class BiopythonOverlayClass:
    child_class = None
    @classmethod
    def cast(cls, entity:bp.Entity.Entity|Self):
        """
        Converts a Bio.PDB object to a bioiain object.
        :param entity: Bio.PDB object to convert.
        :return: bioiain object.
        """

        if isinstance(entity, BiopythonOverlayClass):
            entity.data = deepcopy(entity.data)
            entity.paths = deepcopy(entity.paths)

        entity.__class__ = cls
        entity.base_init()
        if entity.child_class is not None:
            if "child_list" in entity.__dict__.keys():
                entity.__setattr__("child_dict", {})
                for n, child in enumerate(entity.child_list):
                    child.data = deepcopy(entity.data)
                    child.paths = deepcopy(entity.paths)
                    if (isinstance(child, bp.Entity.DisorderedEntityWrapper) and
                            hasattr(entity, "disordered_child_class")):
                        e = entity.disordered_child_class.cast(child)
                    else:
                        e = entity.child_class.cast(child)
                    entity.child_list[n] = e
                    entity.child_dict[n] = e

        if hasattr(entity, "_init"):
            entity._init()
        return entity



    def base_init(self):
        if not hasattr(self, "exporting"):
            self.exporting = ["data", "paths"]
        if not hasattr(self, "data"):
            self.data =  {"info":{
                "name": "_".join([str(i) for i in self.get_full_id()]),
                "o_name": "_".join([str(i) for i in self.get_full_id()]),
                "cls": "",
                "id": 0,
                "repr": "",
            }}
        self.data["info"]["cls"] = self.__class__.__name__
        self.data["info"]["id"] = id(self)
        self.data["info"]["repr"] = repr(self)
        os.makedirs(os.path.abspath("./exports"), exist_ok=True)
        if not hasattr(self, "paths"):
            self.paths = {"export_folder": os.path.abspath("./exports/{}".format(self.data["info"]["o_name"])),
                          "self":None}


    def pass_down(self):
        if "child_list" in self.__dict__.keys():
            for n, child in enumerate(self.child_list):
                child.data = deepcopy(self.data)|child.data
                child.data["info"]["name"] = self.data["info"]["name"] + "_" + str(child.id)
                child.paths = deepcopy(self.paths)|child.paths
                child.paths["self"] = None
                child.paths["parent"] = self.paths["self"]
                child.pass_down()

    def copy_all(self):
        c = self.copy()
        c.data = deepcopy(c.data)
        c.paths = deepcopy(c.paths)
        return c



    def export(self, folder:str|None=None, filename:str|None=None, data:bool=True,
               structure:bool=True, structure_format:str="pdb") -> list[str|None]|str:

        if filename is None:
            filename = self.data["info"]["name"]
        if folder is None:
            folder = self.paths["export_folder"]
        structure_format = structure_format.lower()
        assert structure_format in ["pdb", "cif"]

        paths = []
        os.makedirs(folder, exist_ok=True)

        if structure:
            paths.append(self._export_structure(folder, filename, structure_format))
        if data:
            paths.append(self._export_data(folder, filename))

        if len(paths) == 1:
            return paths[0]
        return paths

    def _export_structure(self, folder, filename, extension="pdb") -> str|None:
        filename = "{}.structure.{}".format(filename, extension)
        filepath = os.path.join(folder, filename)
        self.paths["self"] = filepath
        if extension == "pdb":
            exp = bp.PDBIO()
            exp.set_structure(self)
            exp.save(filepath)
            return filepath
        elif extension == "cif":
            exp = bp.mmcifio.MMCIFIO().set_structure(self)
            exp.save(filepath)
            return filepath
        return None


    def _export_data(self, folder, filename) -> str:

        filename = "{}.data.json".format(filename)
        filepath = os.path.join(folder, filename)
        exp = {e: self.__getattribute__(e) for e in self.exporting if hasattr(self, e)}
        self.paths["data"] = filepath
        with open(filepath, "w") as f:
            f.write(json.dumps(exp, indent=4 ))

        return filepath

    def recover(self):
        if os.path.exists(self.paths["export_folder"]):
            data_path = os.path.join(self.paths["export_folder"], "{}.data.json".format(self.data["info"]["name"]))
            if os.path.exists(data_path):
                log("debug", "recovering data for: {}".format(self.data["info"]["name"]))
                old_data = json.load(open(data_path))
                self.data = self.data | old_data["data"]
                self.paths = self.paths | old_data["paths"]
            else:
                log("debug", "no previous data found for: {}".format(self.data["info"]["name"]))