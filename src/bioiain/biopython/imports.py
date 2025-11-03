import Bio
import Bio.PDB as bp


from .structure import Structure
import Bio.PDB.PDBParser

def loadPDB(file_path:str, name:str=None, quiet=True):
    if name is None:
        name = file_path.split("/")[-1].split(".")[0]
    ext = file_path.split(".")[0]
    if "pdb" in file_path:
        parsed = bp.PDBParser(QUIET=quiet).get_structure(name, file_path)
    elif "cif" in file_path:
        parsed = bp.MMCIFParser(QUIET=quiet).get_structure(name, file_path)
    assert isinstance(parsed, bp.Structure.Structure)
    structure = Structure.cast(parsed)
    return structure


