import Bio.PDB as bp
import os, json

import requests

from ..utilities.strings import clean_string, string_to_list
from ..utilities.logging import log
from .structure import Structure


def loadPDB(file_path:str, name:str=None, quiet=True) -> Structure|None:
    """
    Loads a PDB file into a Structure object.
    :param file_path: Path to the PDB file
    :param name: ID assigned to the structure
    :param quiet: Passed to parser
    :return: Structure object (bioiain)
    """
    if name is None:
        name = file_path.split("/")[-1].split(".")[0]
    ext = file_path.split(".")[-1]
    if "pdb" in ext:
        parsed = bp.PDBParser(QUIET=quiet).get_structure(name, file_path)
    elif "cif" in ext:
        parsed = bp.MMCIFParser(QUIET=quiet).get_structure(name, file_path)
    else:
        log("error", "File format not recognized: {}".format(file_path))
        return None
    assert isinstance(parsed, bp.Structure.Structure)
    structure = Structure.cast(parsed)
    structure.paths["original"] = os.path.abspath(file_path)
    structure.paths["self"] = os.path.abspath(file_path)
    structure.data["info"]["name"] = name
    structure.data["info"]["o_name"] = name
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

    if len(pdb_list) <= 10:
        log("debug", "Codes:", pdb_list)
    else:
        log("debug", "Codes:", len(pdb_list))

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




def recover(name, export_folder="./exports", download_dir="./data", download=True):
    log(2, "Recovering structure: {}".format(name))
    pdb_code = name.split("_")[0].upper()
    log("debug", "PDB code: {}".format(pdb_code))

    try:
        exported_folders = os.listdir(export_folder)
    except:
        log("warning", "Export folder does not exist: {}".format(os.path.abspath(export_folder)))
        exported_folders = None

    if exported_folders != None:
        if pdb_code in exported_folders:
            exported_folder = os.path.join(export_folder, pdb_code)
            json_path  = os.path.join(exported_folder, name+".data.json")
            print(os.path.abspath(json_path))
            with open(json_path, "r") as f:
                data = json.load(f)
            #print(data)

            structure = loadPDB(data["paths"]["original"])
            for k, v in data.items():
                setattr(structure, k, v)
            #print(structure.data)
            #print(structure.paths)
            return structure



        else:
            log("warning", "Export folder of {} not found".format(pdb_code))

    return None


def read_mmcif(file_path, subset:list|str=None, exclude:list|str=None) -> dict:
    from ..utilities.strings import str_to_list_with_literals
    data = {}
    with open(file_path, "r") as f:
        n = -1
        in_loop = False
        group_key = None
        loop_keys = None
        loop_values = None
        next_line = None
        eof = False
        multi_line = False
        multi_cached = None
        multi_delimiters = ["\'","\"",";"]
        while not eof:
            n+=1
            line = next_line
            next_line = next(f, None)
            #print(f"{n}>>{repr(line)}<<")



            if line is None:
                continue
            if next_line is None:
                eof = True
            if line.startswith("#"):
                if in_loop:
                    print("LOOP END")
                else:
                    print("SECTION END")
                in_loop = False
                group_key = None
                continue

            if line.startswith("loop_"):
                in_loop = True
                loop_keys = []
                loop_values = []
                group_key = []
                print("LOOP START")
                continue

            if multi_line:
                #print(multi_line)
                if multi_cached is None:
                    multi_cached = ""
                if line == "\n":
                    multi_cached += line
                    continue
                # print(repr(line))
                # print([l for l in line])
                # if line.startswith(";\n"):
                #     break
                if line[-2] in multi_delimiters:
                    #print(repr(line))
                    #print(repr(line[-3:-3]))
                    multi_line = False
                    if in_loop:
                        v.append(multi_cached)
                        print("MULTI LINE END")
                    else:
                        print("MULTI LINE END (in-loop)")
                    multi_cached = None

                else:
                    if line[0] in multi_delimiters:
                        multi_cached += line[1:]
                    else:
                        multi_cached += line



            if not in_loop:

                if not multi_line:
                    #print(line)
                    line_list = str_to_list_with_literals(line)
                    #print(repr(line_list))
                    if len(line_list) == 0:
                        continue
                    group_key = line_list[0].split(".")[0]
                    try:
                        k = line_list[0].split(".")[1]
                    except IndexError:

                        k = group_key
                        group_key = None

                    if len(line_list) == 1:
                        v = []
                    else:
                        v = line_list[1:]

                if next_line[0] in multi_delimiters:
                    print("MULTI LINE START")
                    multi_line = True
                    continue

                if group_key is None:
                    print(repr(line))
                    print(line_list)
                    log("warning", "no key-value structure found in:", repr(line))
                if not multi_line:
                    if group_key not in data.keys():
                        data[group_key] = {}
                    data[group_key][k] = " ".join(v)


            else: # IN LOOP
                if line.startswith("_"):
                    group_key.append(line.split(".")[0])
                    loop_keys.append(line.split(" ")[0].split(".")[1])
                    continue
                else:
                    if next_line[0] in multi_delimiters:
                        print("MULTI LINE START (in-loop)")
                        multi_line = True
                        continue
                    loop_values.append(line)

    json.dump(data, open("out.json", "w"), indent=4)



