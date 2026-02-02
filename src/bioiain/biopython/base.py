import builtins
import os
from copy import deepcopy

import Bio.PDB as bp
import json


from ..utilities.logging import log
from typing_extensions import Self, LiteralString


class BiopythonOverlayClass:
    child_class = None

    def _getitem(self, item):
        l = [c for c in self.child_dict.values() if c.id == item]
        if len(l) == 0:
            raise KeyError(item)
        elif len(l) == 1:
            return l[0]
        else:
            raise Exception("Multiple children with matching id")



    @classmethod
    def cast(cls, entity:bp.Entity.Entity|Self, *args, **kwargs) -> bp.Entity.Entity:
        """
        Converts a Bio.PDB object to a bioiain object.
        :param entity: Bio.PDB object to convert.
        :return: bioiain object.
        """

        bio_id = None
        if isinstance(entity, BiopythonOverlayClass):
            entity.data = deepcopy(entity.data)
            entity.paths = deepcopy(entity.paths)
        else:
            bio_id = entity.get_id()

        entity.__class__ = cls
        entity.base_init()

        if bio_id is not None:
            entity.data["info"]["bio_id"] = bio_id

        if entity.child_class is not None:
            if "child_list" in entity.__dict__.keys():
                entity.__setattr__("child_dict", {})
                for n, child in enumerate(entity.child_list):
                    for attr in entity.exporting:
                        setattr(child, attr, deepcopy(getattr(entity, attr)))

                    if (isinstance(child, bp.Entity.DisorderedEntityWrapper) and
                            hasattr(entity, "disordered_child_class")):
                        e = entity.disordered_child_class.cast(child)
                    else:
                        e = entity.child_class.cast(child)
                    entity.child_list[n] = e
                    entity.child_dict[n] = e

        if hasattr(entity, "_init"):
            entity._init(*args, **kwargs)

        return entity



    def base_init(self):
        if not hasattr(self, "exporting"):
            self.exporting = ["data", "paths", "flags"]
        if not hasattr(self, "data"):
            self.data =  {"info":{
                "name": "_".join([str(i) for i in self.get_full_id()]),
                "o_name": "_".join([str(i) for i in self.get_full_id()]),
                "cls": "",
                "id": 0,
                "repr": "",
                "bio_id": None,
            }}
        self.data["info"]["cls"] = self.__class__.__name__
        self.data["info"]["id"] = id(self)
        self.data["info"]["repr"] = repr(self)
        os.makedirs(os.path.abspath("./exports"), exist_ok=True)
        if not hasattr(self, "paths"):
            self.paths = {"export_folder": os.path.abspath("./exports/{}".format(self.data["info"]["o_name"])),
                          "self":None}
        if not hasattr(self, "flags"):
            self.flags = {}

        self.flags["base_init"] = True
        self._atoms = None


    def pass_down(self):
        if "child_list" in self.__dict__.keys():
            for n, child in enumerate(self.child_list):
                for attr in self.exporting:
                    setattr(child, attr, deepcopy(getattr(self, attr))|getattr(child, attr))

                child.data["info"]["name"] = self.data["info"]["name"] + "_" + str(child.id)
                child.paths["self"] = None
                child.paths["parent"] = self.paths["self"]
                child.pass_down()

    def copy_all(self):
        c = self.copy()
        for attr in self.exporting:
            setattr(c, attr, deepcopy(getattr(c, attr)))
        return c

    def get_name(self):
        return self.data["info"]["name"]

    def has_flag(self, key, value=None):
        if value is None:
            return key in self.flags
        else:
            if key in self.flags:
                if value == self.flags[key]:
                    return True
        return False

    def add_flag(self, key, value=None, export=True):
        self.flags[key] = value
        if export:
            self.export()

    def export(self, folder:str|None=None, filename:str|None=None, data:bool=True,
               structure:bool=True, structure_format:str="cif") -> list[str|None]|str:

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

        self.flags["exported"] = True

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
            exp = bp.mmcifio.MMCIFIO()
            exp.set_structure(self)
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

    @classmethod
    def recover(cls, *args, data_path:str|LiteralString=None, load_structure:bool=True, temp_id="recovering", **kwargs):
        if not data_path.endswith(".data.json"):
            data_path += ".data.json"
        try:
            self = cls()
        except:
            try:
                self = cls(*args, **kwargs)
            except:
                self = cls(temp_id, *args, *kwargs)

        self.base_init()


        if load_structure:
            from .imports import loadPDB
            data = self.read_data_file(data_path)["data"]
            struc_path = data_path.replace(".data.json", ".structure.cif")
            struc = loadPDB(struc_path, data["info"]["o_name"])
            child = struc
            deepness = 0
            while not issubclass(self.__class__, child.__class__):
                try:
                    child = child[0]
                    deepness += 1
                except:
                    break
            if issubclass(self.__class__, child.__class__):

                self = self.cast(child)
            else:
                raise Exception("Failed to load structure")

        #print(self)
        r = self._recover(data_path=data_path)
        if r is None:
            log("warning", f"Failed recovery for: {data_path}")
            return None
        #print(self.id)
        self.parent = {"child_dict": {}}
        #self.__setattr__("_id", str(self.get_name()))
        # for path in self.paths:
        #     path.replace(temp_id, self.get_name())
        #self.pass_down()
        #self.__setattr__("id", self.get_name())

        #print(self)
        return self


    def _recover(self, data_path=None):
        if data_path is None:
            name = self.data["info"]["name"]
            data_path = os.path.join(self.paths["export_folder"], f"{name}.data.json")
        else:
            name = data_path.split(".")[0]
        if data_path is None:
            log("warning", "No data path provided for recover()")
            return None

        elif os.path.exists(data_path):
            log("debug", f"recovering data for: {name}")
            #print(data_path)
            old_data = json.load(open(data_path))
            #print(self.data)
            #print(json.dumps(old_data["data"], indent=4))
            for attr in self.exporting:
                self.__getattribute__(attr).update(old_data[attr])

            #print(json.dumps(self.data, indent=4))
            self.add_flag("recovered", True)
            self.export()
            return self
        else:
            log("warning", f"No previous data found for: {name}")
            return None

    def read_data_file(self, filepath):
        return json.load(open(filepath))


class Entity(bp.Entity.Entity, BiopythonOverlayClass):
    child_class = None

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)


class BioiainObject(object):
    pass
