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



class MMCIF(object):
    def __init__(self, data):
        self.data = data

    def save(self, path):
        json.dump(self.data, open(path, "w"), indent=4)

    def __getitem__(self, key):
        index = None
        subkey = None
        if type(key) in [list, tuple]:
            if len(key) == 1:
                key = key[0]
            if len(key) == 2:
                key, index = key
            elif len(key) == 3:
                key, index, subkey = key

        if index is None:
            index = 0


        ret = self.data[key][index]

        if subkey is not None:
            ret = ret[subkey]

        return ret









def read_mmcif(file_path, output_folder="headers", subset:list|str=None, exclude:list|str=None) -> MMCIF:
    from ..utilities.strings import str_to_list_with_literals
    data = {}
    name = os.path.basename(file_path).split(".")[0]
    os.makedirs(output_folder, exist_ok=True)
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
        multi_delimiter = None
        multi_delimiters = [";"]
        looping = False
        while not eof:
            n+=1
            line = next_line
            next_line = next(f, None)
            if line is None:
                continue
            if next_line is None:
                eof = True
            if line.startswith("#"):
                in_loop = False
                looping = False
                group_key = None
                multi_cached = None
                continue

            if line.startswith("loop_"):
                in_loop = True
                loop_keys = []
                loop_values = []
                group_key = []
                continue

            if multi_line:
                if multi_cached is None:
                    multi_cached = ""
                if line.replace("\n", "").strip() == "":
                    multi_cached += line
                    continue
                if type(multi_cached) is list:
                    line_list, open_lit = str_to_list_with_literals(line, check_open_literal=True )

                if line.replace("\n", "").strip()[-1] in multi_delimiters:
                    line = line.replace("\n", "").strip()[:-2]
                    multi_line = False
                    multi_delimiter = None
                    if not in_loop:
                        v.append(multi_cached)
                        multi_cached = None
                        continue
                    else:
                        loop_values[-1].append(multi_cached)


                else:
                    if type(multi_cached) is str:
                        if line[0] in multi_delimiters:
                            multi_cached += line[1:]
                        else:
                            multi_cached += line
                    elif type(multi_cached) is list:
                        if open_lit:
                            multi_cached[-1] += line_list[0]
                            multi_cached.extend(line_list[1:])

                        else:
                            multi_cached.extend(line_list)
                    else:
                        raise TypeError("Bioiain mmcif Parser error")



            if not in_loop:

                if not multi_line:
                    line_list, open_lit = str_to_list_with_literals(line, check_open_literal=True)
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
                if open_lit:
                    multi_line = True
                    multi_cache = line_list[-1]
                    multi_delimiter = line_list[-1][0]
                if next_line[0] in multi_delimiters and not multi_line:
                    multi_line = True
                    multi_delimiter = next_line[0]
                    continue

                if group_key is None:
                    #print("multi line", multi_line)
                    #print("line",repr(line))
                    #print("line_list",line_list)
                    #print(n)
                    log("warning", f"No key-value structure found in line {n}:", repr(line))
                if not multi_line:
                    if exclude is not None:
                        if group_key in exclude:
                            continue
                    if subset is not None:
                        if group_key not in subset:
                            continue
                    if group_key not in data.keys():
                        data[group_key] = [{}]
                    #print(f"{group_key}.{k} ->", " ".join(v))
                    data[group_key][0][k] = " ".join(v)


            else: # IN LOOP
                if not multi_line:
                    if line.startswith("_") and not looping:
                        group_key.append(line.split(".")[0])
                        loop_keys.append(line.split(" ")[0].split(".")[1])
                        continue
                    else:
                        looping = True
                    if multi_cached is None:
                        if len(loop_values) > 0:
                            if len(loop_values[-1]) < len(loop_keys):
                                loop_values[-1].extend(str_to_list_with_literals(line))
                            else:
                                loop_values.append(str_to_list_with_literals(line))
                        else:
                            loop_values.append(str_to_list_with_literals(line))
                    else:
                        multi_cached = None
                if (next_line[0] in multi_delimiters) and (not multi_line):
                    multi_line = True
                    multi_delimiter = next_line[0]
                    multi_cache = []
                    continue

                if next_line[0] == "#":
                    try:
                        assert len(set(group_key)) == 1
                    except AssertionError:
                        log("warning", f"Multiple keys structure found in line {n}:", repr(line))
                        continue
                    group_key = group_key[0]
                    if exclude is not None:
                        if group_key in exclude:
                            continue
                    if subset is not None:
                        if group_key not in subset:
                            continue
                    if group_key not in data.keys():
                        data[group_key] = []
                    for i, l in enumerate(loop_values):
                        try:
                            assert len(l) == len(loop_keys)
                        except AssertionError:
                            log("warning", f"Missmatch on key-value numbers in group {group_key}, element {i}")
                            continue
                        d = {k:v for k, v in zip(loop_keys, l)}
                        data[group_key].append(d)
                    #print(f"{group_key} ->", f"list of length: {len(loop_values)}")
    save_path = os.path.join(output_folder, f"{name}.header.json")
    log("debug","Headers saved to:", os.path.abspath(save_path))
    mmcif = MMCIF(data)
    mmcif.save(save_path)
    return mmcif



