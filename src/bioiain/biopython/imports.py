import Bio.PDB as bp
import os

import requests

from ..utilities.strings import clean_string, string_to_list
from ..utilities.logging import log
from .structure import Structure


def loadPDB(file_path:str, name:str=None, quiet=True) -> Structure:
    """
    Loads a PDB file into a Structure object.
    :param file_path: Path to the PDB file
    :param name: ID assigned to the structure
    :param quiet: Passed to parser
    :return: Structure object (bioiain)
    """
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
    structure.paths["original"] = os.path.abspath(file_path)
    return structure


def downloadPDB(data_dir:str, list_name:str, pdb_list:list=None, file_path:str = None, file_format="pdb",
                overwrite:bool=False) -> str:
    """
    Downloads a list of PDB files into a folder of given name within the data_dir. Creates a file containing all
    the file in the data_dir
    :param data_dir: Directory to create the download folder
    :param list_name: Name of the folder to download files to
    :param pdb_list: List of PDB codes to download
    :param file_path: (optional) file with PDB codes, separated by comma or new lines, extends pdb_list
    :param file_format: PDB / CIF, extension of downloaded files
    :param overwrite: True to download existing pdb files and overwrite
    :return: Path to folder containing downloaded files
    """
    log("debug", "Downloading PDB files...")
    file_format = file_format.lower()
    assert file_format in ["pdb", "cif"]
    if pdb_list is None:
        pdb_list = []
    if file_path is not None:
        with open(file_path) as f:
            for line in f:
                new = string_to_list(line, delimiter=",")
                for n in new:
                    n = clean_string(n)
                    pdb_list.append(n)
    pdb_list = sorted(list(set([p.upper() for p in pdb_list])))

    log("debug", "Codes:", pdb_list)

    os.makedirs(data_dir, exist_ok=True)
    list_folder = os.path.join(data_dir, list_name)
    os.makedirs(list_folder, exist_ok=True)
    link_file = "{}_({}).txt".format(list_name, file_format)
    with open(os.path.join(data_dir, link_file) , "w") as f:
        for pdb in pdb_list:
            if file_format == "pdb":
                f.write("https://files.rcsb.org/download/{}.pdb\n".format(pdb))
            elif file_format == "cif":
                f.write("https://files.rcsb.org/download/{}.cif\n".format(pdb))
    log("debug", "Generated links at: {}".format(os.path.join(data_dir, link_file) ))

    with open(os.path.join(data_dir, link_file)) as f:
        counter = 0
        failed_counter = 0
        skipped_counter = 0
        for line in f:
            line = line.replace("\n", "")
            f_name = line.split("/")[-1]
            if os.path.exists(os.path.join(list_folder, f_name)) and not overwrite:
                skipped_counter += 1
                continue
            url = line
            log("debug", "...Downloading {}".format(url), end="\r")
            response = requests.get(url)
            if response.status_code != 200:
                log("Error", "Failed to download from:", line)
                failed_counter += 1
            else:
                with open(os.path.join(list_folder, f_name), "w") as f:
                    f.write(response.text)
                counter += 1
    log("debug", "{} files downloaded, {} failed, {} skipped".format(counter, failed_counter, skipped_counter))
    return list_folder





