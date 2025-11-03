import Bio.PDB as bp
from ..utilities.logging import log
from .structure import Structure


def loadPDB(file_path:str, name:str=None, quiet=True):
    if name is None:
        name = file_path.split("/")[-1].split(".")[0]
    ext = file_path.split(".")[0]
    if "pdb" in file_path:
        parsed = bp.PDBParser(QUIET=quiet).get_structure(name, file_path)
    elif "cif" in file_path:
        parsed = bp.MMCIFParser(QUIET=quiet).get_structure(name, file_path)
    else:
        log("error", "File format not recognized: {}".format(file_path))
        return None
    assert isinstance(parsed, bp.Structure.Structure)
    structure = Structure.cast(parsed)
    return structure


