import os, json

from ..utilities import *



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
                        log("warning", f"Multiple keys structure found in line {n}:", repr(line))
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

def write_atoms(atoms, file_path, name=None, include_unused=True, include_misc=True, preserve_ids=False,
                mode="w", key="_atom_site") -> str:

    labels = atoms[0]._mmcif_dict(include_unused=include_unused,  include_misc=include_misc).keys()

    if not file_path.endswith(".cif"):
        file_path += ".cif"
    os.makedirs(os.path.dirname(file_path), exist_ok=True)
    try:
        with open(file_path, mode) as f:
            if mode == "w" and name is not None:
                f.write(f"data_{name}\n")
            f.write("#\n")
            f.write("loop_\n")

            for l in labels:
                f.write(f"{key}.{l}\n")

            for n, a in enumerate(atoms):
                d = a._mmcif_dict(include_unused=include_unused, include_misc=include_misc)
                if not preserve_ids:
                    d["id"] = f"{n:4d}"
                f.write("  ".join(d.values()) + "\n")

        return file_path
    except Exception as e:
        log("error", f"Atom write interrupted, deleting corrupted file: {file_path}")
        os.remove(file_path)
        log("error", "File deleted successfully!")
        raise e


def write_dict(data, label, file_path, name=None, mode="w"):

    if not file_path.endswith(".cif"):
        file_path += ".cif"
    if not label.startswith("_"):
        label = "_" + label
    try:
        with open(file_path, mode) as f:
            if mode == "w" and name is not None:
                f.write(f"data_{name}\n")
            f.write("#\n")

            for k, v in data.items():
                f.write(f"{label}.{k}   {str(v)}\n")
        return file_path
    except Exception as e:
        log("error", f"Dict writing to mmcif failed, deleting corrupted file: {file_path}")
        os.remove(file_path)
        raise e

def write_list(data, label, file_path, name=None, mode="w"):
    pass


