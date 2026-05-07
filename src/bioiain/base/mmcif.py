import os, json, requests

from ..utilities import *

from ..utilities.strings import *

def downloadPDB(data_dir:str, list_name:str, pdb_list:list=None, file_path:str = None, file_format="cif",
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
    log("header", "Downloading PDB files...")
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
        log(1, "N codes:", pdb_list)
    else:
        log(1, "N codes:", len(pdb_list))

    os.makedirs(data_dir, exist_ok=True)
    list_folder = os.path.join(data_dir, list_name)
    os.makedirs(list_folder, exist_ok=True)
    link_file = "{}_{}.link.list".format(list_name, file_format)
    with open(os.path.join(data_dir, link_file) , "w") as f:
        for pdb in pdb_list:
            if file_format == "pdb":
                f.write("https://files.rcsb.org/download/{}.pdb\n".format(pdb))
            elif file_format == "cif":
                f.write("https://files.rcsb.org/download/{}.cif\n".format(pdb))
    # log("debug", "Generated links at: {}".format(os.path.join(data_dir, link_file) ))

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
    log(1, "{} files downloaded, {} failed, {} skipped".format(counter, failed_counter, skipped_counter))
    return list_folder

class MMCIF(object):
    def __init__(self, data, cif_path=None):
        self.data = data
        self.cif_path = cif_path

    def save(self, path):
        json.dump(self.data, open(path, "w"), indent=4)

    def __getitem__(self, key):
        index = None
        subkey = None
        key = key.split(".")
        #print(key)

        if type(key) is str:
            pass
        elif len(key) == 1:
            key = key[0]
        elif len(key) == 2:
            key, subkey = key

        elif len(key) == 3:
            key, index, subkey = key
        else:
            log("warning", "Invalid key {}".format(key))

        if not key.startswith("_"):
            key = "_" + key
        try:
            d = self.data[key]
            if len(d) == 1 and index is None:
                index = 0
            #print(key, index, subkey)
            if index is None and subkey is None:
                ret = [v for v in d]
            elif index is None:
                #print([v.keys() for n, v in enumerate(d)])
                ret = [v[subkey] for v in d]
            elif subkey is None:
                ret = d[index]
            else:
                ret = d[index][subkey]
        except KeyError as e:
            #print(ret)
            log("warning", f"Key not found: {e} ({key}.{index}.{subkey}) in {self.cif_path}")
            return None

        return ret

    def __call__(self, *args):
        entry = ".".join([*args])
        return self[entry]





def read_mmcif(file_path, output_folder=None, subset:list|str=None, exclude:list|str=None) -> MMCIF:
    from ..utilities.strings import str_to_list_with_literals
    data = {}
    name = os.path.basename(file_path).split(".")[0]
    if output_folder is not None:
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
                next_line = "#"
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
                #print("LOOP START")
                continue
            #print("\nLINE:")
            #print(repr(line))
            #print("group_key:", repr(group_key))
            #print("multi:", multi_line)
            #print("loop:", in_loop)
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
                    open_lit = False
                    multi_delimiter = None
                    #print("MULTI LINE END")
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
                    group_key = line_list[0].split(".")[0].replace("\n", "").strip()
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
                    #print(multi_cached)
                    #print(line_list)
                    #print("MULTI_LINE START (OPEN-LIT)")
                    multi_line = True
                    multi_cache = line_list[-1]
                    multi_delimiter = line_list[-1][0]

                if next_line[0] in multi_delimiters and not multi_line:
                    #print("MULTI_LINE START")
                    multi_line = True
                    multi_delimiter = next_line[0]
                    continue

                if group_key is None:
                    #print("multi line", multi_line)
                    #print("line",repr(line))
                    #print("line_list",line_list)
                    #print(n)
                    if n != 1:
                        log("warning", f"No key-value structure found in line {n}:", repr(line), f"\n  (In file: {file_path})")

                    else:
                        #log("debug", "Parsing:", line.replace("\n", "").strip())
                        pass

                if not multi_line and group_key is not None:
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
                #print("MULTI:", multi_line)
                if not multi_line:
                    if line.startswith("_") and not looping:
                        group_key.append(line.split(".")[0].replace("\n", "").strip())
                        loop_keys.append(line.split(" ")[0].split(".")[1].replace("\n", "").strip())
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
                        print(group_key)
                        log("error", f"Multiple keys structure found in line {n}:", repr(line))
                        exit()
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

    mmcif = MMCIF(data, cif_path=file_path)
    if output_folder is not None:
        save_path = os.path.join(output_folder, f"{name}.header.json")
        log("debug","Headers saved to:", os.path.abspath(save_path))
        mmcif.save(save_path)
    return mmcif


def write_pdb_atoms(atoms, file_path, mode="w", end=True):

    if not file_path.endswith(".pdb"):
        file_path += ".pdb"
    os.makedirs(os.path.dirname(file_path), exist_ok=True)

    try:
        with open(file_path, mode) as f:

            rn = 1
            for an, atom in enumerate(atoms):
                if (atom.chain != atoms[an-1].chain and an != 0):
                    rn += 1
                    ter = atoms[an-1].pdb_string(rn)
                    ter = "TER   "+ter[6:27]
                    f.write(ter + "\n")
                s = atom.pdb_string(rn)
                f.write(s + "\n")
                rn += 1
                if an == len(atoms) - 1:
                    rn += 1
                    ter = atoms[an - 1].pdb_string(rn)
                    ter = "TER   " + ter[6:27]
                    f.write(ter + "\n")
            if end:
                f.write("END\n")


    except Exception as e:
        log("error", f"Atom write (pdb) interrupted, deleting corrupted file: {file_path}")
        os.remove(file_path)
        log("error", "File deleted successfully!")
        raise e




def write_atoms(atoms, file_path, name=None, include_misc=True, preserve_ids=False,
                mode="w", key="_atom_site") -> str:
    if len(atoms) == 0:
        return None
    labels = atoms[0]._mmcif_dict( include_misc=include_misc).keys()

    if not file_path.endswith(".cif"):
        file_path += ".cif"
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    if name is None:
        name = os.path.basename(file_path).split(".")[0]
    try:
        with open(file_path, mode) as f:
            if mode == "w":
                f.write(f"data_{name}\n")
            f.write("#\n")
            f.write("loop_\n")

            for l in labels:
                f.write(f"{key}.{l}\n")
            n = 1
            for a in atoms:
                d = a._mmcif_dict(include_misc=include_misc)
                if not preserve_ids:
                    d["id"] = f"{n:4d}"
                f.write("  ".join(d.values()) + "\n")
                n+=1

        return file_path
    except Exception as e:
        log("error", f"Atom write interrupted, deleting corrupted file: {file_path}")
        os.remove(file_path)
        log("error", "File deleted successfully!")
        raise e


def write_dict(data, label, file_path, name=None, mode="w"):

    if len(data) == 0:
        return file_path
    if not file_path.endswith(".cif"):
        file_path += ".cif"
    if not label.startswith("_"):
        label = "_" + label
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    try:
        with open(file_path, mode) as f:
            if mode == "w" and name is not None:
                f.write(f"data_{name}\n")
            f.write("#\n")

            for k, v in data.items():
                if v is None:
                    f.write(f"{label}.{k}   ?\n")
                elif f"{label}.{k}" in quoted_headers or should_be_quoted(v):
                    f.write(f"{label}.{k}   '{str(v)}'\n")
                else:
                    f.write(f"{label}.{k}   {str(v)}\n")
        return file_path
    except Exception as e:
        log("error", f"Dict writing to mmcif failed, deleting corrupted file: {file_path}")
        os.remove(file_path)
        raise e

def write_dict_list(data, label, file_path, name=None, mode="w", **kwargs):

    if len(data) == 0:
        return file_path
    if not file_path.endswith(".cif"):
        file_path += ".cif"
    if not label.startswith("_"):
        label = "_" + label

    log(3, "Writing dict list to:", file_path)
    log(4, "Label:", label)

    keys = ["n"]

    if type(data) is dict:
        data = data.values()
    else:
        assert type(data) in (list, tuple)

    get_dict=False

    if type(data[0]) is not dict:
        get_dict=True

    for d in data:
        if d is not None:
            if get_dict:
                keys.extend(d._mmcif_dict(**kwargs).keys())
            else:
                keys.extend(d.keys())
            break


    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    try:
        with open(file_path, mode) as f:
            if mode == "w" and name is not None:
                f.write(f"data_{name}\n")
            f.write("#\n")
            f.write("loop_\n")

            for k in keys:
                f.write(f"{label}.{k}\n")

            for n, d in enumerate([d for d in data if d is not None]):
                if get_dict:
                    d = d._mmcif_dict(**kwargs)
                    f.write(f"{n:4d}  "+"  ".join([quote_if_necessary(v) for v in d.values()]) + "\n")

        return file_path
    except Exception as e:
        log("error", f"Dict writing to mmcif failed, deleting corrupted file: {file_path}")
        os.remove(file_path)
        raise e


def quote_if_necessary(value):
    if should_be_quoted(value):
        return f"'{value}'"
    else:
        return str(value)


def should_be_quoted(value):
    if value is None:
        return False
    value = str(value)
    if value.strip() == "":
        return True
    if value.startswith("'") and value.endswith("'"):
        return False
    if "'" in value:
        return True
    if value.startswith("<") or value.endswith(">"):
        return True

    return False



quoted_headers = [
    "_symmetry.space_group_name_H-M",
]
