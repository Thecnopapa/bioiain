import builtins
import os
import Bio.PDB as bp
import json
from ..utilities.logging import log



class BiopythonOverlayClass:
    child_class = None
    @classmethod
    def cast(cls, entity:bp.Entity.Entity):
        """
        Converts a Bio.PDB object to a bioiain object.
        :param entity: Bio.PDB object to convert.
        :return: bioiain object.
        """
        entity.__class__ = cls
        if entity.child_class is not None:
            if "child_list" in entity.__dict__.keys():
                entity.__setattr__("child_dict", {})
                for n, child in enumerate(entity.child_list):
                    if (isinstance(child, bp.Entity.DisorderedEntityWrapper) and
                            hasattr(entity, "disordered_child_class")):
                        e = entity.disordered_child_class.cast(child)
                    else:
                        e = entity.child_class.cast(child)
                    entity.child_list[n] = e
                    entity.child_dict[n] = e
        entity.base_init()
        if hasattr(entity, "init"):
            entity.init()
        return entity




    def base_init(self):
        if not hasattr(self, "exporting"):
            self.exporting = ["data", "paths"]
        if not hasattr(self, "data"):
            self.data =  {"info":{"name":"_".join([str(i) for i in self.get_full_id()])}}
        if not hasattr(self, "paths"):
            self.paths = {"export_folder": os.path.abspath("./exports")}
            os.makedirs(self.paths["export_folder"], exist_ok=True)



    def export(self, folder:str|None=None, filename:str|None=None, data:bool=True,
               structure:bool=True, structure_format:str="pdb") -> list[str|None]|str:

        if filename is None:
            filename = self.data["info"]["name"]
        if folder is None:
            folder = os.path.join(self.paths["export_folder"], filename.split("_")[0])
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
        filename = "{}.{}".format(filename, extension)
        filepath = os.path.join(folder, filename)
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
        with open(filepath, "w") as f:
            f.write(json.dumps(exp, indent=4 ))
        return filepath
