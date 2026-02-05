import os, sys, subprocess, json
from ..utilities.logging import log
from ..utilities.sequences import d3to1

def ss_to_index(ss):

    if ss == "H":
        return 0
    if ss == "B":
        return 1
    if ss == "G":
        return 2
    if ss == "I":
        return 3
    if ss == "P":
        return 4
    if ss == "E":
        return 5
    if ss == "T":
        return 6
    if ss == "S":
        return 7
    if ss == "-" or ss == " " or ss=="#":
        return 8
    return 8



def index_to_ss(ss):

    if ss == 0:
        return "H"
    if ss == 1:
        return "B"
    if ss == 2:
        return "G"
    if ss == 3:
        return "I"
    if ss == 4:
        return "P"
    if ss == 5:
        return "E"
    if ss == 6:
        return "T"
    if ss == 7:
        return "S"
    if ss == 8:
        return "#"
    return "-"



class DSSP:
    def __init__(self, dssp_cmd="dssp", folder="/tmp/bioiain/dssp", force=False):

        self.dssp_cmd = dssp_cmd
        self.folder = folder
        os.makedirs(self.folder, exist_ok=True)


    def eval(self, file_path, name=None, force=False):
        if name is None:
            name = os.path.basename(file_path).split(".")[0]
        fname = os.path.join(file_path, name+".dssp")
        out_path = os.path.join(self.folder, name+".dssp")
        if os.path.exists(out_path) and not force:
            return fname

        out_path = self._run_dssp(file_path, out_path)
        #dssp_dict = self._read_dssp(out_path)


    def _run_dssp(self, file_path, out_path, force=False):
        cmd = [self.dssp_cmd, "--output-format", "mmcif", "-v",  f"{file_path}", f"{out_path}"]
        print("$", " ".join(cmd))
        subprocess.run(cmd)
        return out_path


    def _read_dssp(self, out_path):

        real_chains = list(label_dict.keys())
        halucination_pointer = {}
        # print("REAL", real_chains)
        dssp_dict = {}
        with open(f"{raw_folder}/{name}.dssp", "r") as f:
            start_sum = False
            start_bridge = False
            for line in f:
                if "_dssp_struct_summary" in line:
                    start_sum = True
                    continue
                if "_dssp_struct_bridge" in line:
                    start_bridge = True
                    continue

                if "#" in line:
                    start_sum = False
                    if start_bridge:
                        break
                    continue
                if start_sum:
                    # print(line)
                    line = line.split(" ")
                    line = [l for l in line if l != ""]
                    # print(line)
                    res = line[2]
                    dch = line[1]
                    if dch not in dssp_dict.keys():
                        dssp_dict[dch] = {}
                    ch = real_chains[list(dssp_dict.keys()).index(dch)]

                    # print("DSSP:", dch, "REAL:", ch)
                    resn = line[3].upper()
                    ss = line[4].replace(".", "#")

                    if ch not in label_dict.keys():
                        continue
                    try:
                        dssp_dict[dch][len(dssp_dict[dch])] = {"res": res, "resn": d3to1[resn], "resn3": resn,
                                                               abbreviation: ss, "label": ss_to_index(ss)}
                    except:
                        dssp_dict[dch][len(dssp_dict[dch])] = {"res": res, "resn": None, "resn3": resn,
                                                               abbreviation: ss, "label": ss_to_index(ss)}
                elif start_bridge:
                    # print(line)
                    line = line.split(" ")
                    line = [l for l in line if l != ""]
                    # print(line)
                    dch = line[3].upper()
                    rch = line[5].upper()

                    if dch not in halucination_pointer.keys():
                        halucination_pointer[dch] = rch
        # print(halucination_pointer)
        [log(4, "DSSP:", dch, "->>", "REAL:", rch) for dch, rch in halucination_pointer.items()]

        for dch, rch in halucination_pointer.items():
            label_dict[rch] = dssp_dict[dch]

        # print("DSSP", dssp_dict.keys())
        if not start_bridge:
            log("error", "Error reading DSSP file for:", name)


def run_dssp(name, label_dict, data_folder, save_folder, raw_folder, abbreviation, force=False, dssp_command="dssp"):
    os.makedirs(raw_folder, exist_ok=True)
    os.makedirs(save_folder, exist_ok=True)

    if os.path.exists(f"{save_folder}/{name}.labels.json") and not force:
        log(3, "Label already exists")
        return True

    cmd = [dssp_command, "--output-format", "mmcif", "-v", f"{data_folder}/{name}.cif",
           f"{raw_folder}/{name}.dssp"]
    log(3, " ".join(cmd))
    subprocess.run(cmd)

    real_chains = list(label_dict.keys())
    halucination_pointer = {}
    #print("REAL", real_chains)
    dssp_dict = {}
    with open(f"{raw_folder}/{name}.dssp", "r") as f:
        start_sum = False
        start_bridge = False
        for line in f:
            if "_dssp_struct_summary" in line:
                start_sum = True
                continue
            if "_dssp_struct_bridge" in line:
                start_bridge = True
                continue


            if "#" in line:
                start_sum = False
                if start_bridge:
                    break
                continue
            if start_sum:
                #print(line)
                line = line.split(" ")
                line = [l for l in line if l != ""]
                #print(line)
                res = line[2]
                dch = line[1]
                if dch not in dssp_dict.keys():
                    dssp_dict[dch] = {}
                ch = real_chains[list(dssp_dict.keys()).index(dch)]

                #print("DSSP:", dch, "REAL:", ch)
                resn = line[3].upper()
                ss = line[4].replace(".", "#")

                if ch not in label_dict.keys():
                    continue
                try:
                    dssp_dict[dch][len(dssp_dict[dch])] = {"res": res, "resn": d3to1[resn], "resn3":resn, abbreviation: ss, "label": ss_to_index(ss)}
                except:
                    dssp_dict[dch][len(dssp_dict[dch])] = {"res": res, "resn": None, "resn3":resn, abbreviation: ss, "label": ss_to_index(ss)}
            elif start_bridge:
                # print(line)
                line = line.split(" ")
                line = [l for l in line if l != ""]
                # print(line)
                dch = line[3].upper()
                rch = line[5].upper()

                if dch not in halucination_pointer.keys():
                    halucination_pointer[dch] = rch
    #print(halucination_pointer)
    [log(4, "DSSP:", dch, "->>", "REAL:", rch) for dch, rch in halucination_pointer.items()]

    for dch, rch in halucination_pointer.items():
        label_dict[rch] = dssp_dict[dch]




    #print("DSSP", dssp_dict.keys())
    if not start_bridge:
        log("error", "Error reading DSSP file for:", name)

    with open(f"{save_folder}/{name}.labels.json", "w") as f:
        f.write(json.dumps(label_dict, indent=4))

    return True

