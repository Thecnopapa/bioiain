import os, json
from sys import orig_argv

import numpy as np

from .space_groups import dictio_space_groups
from ..biopython import read_mmcif
from ..utilities.logging import log
from ..utilities.strings import string_to_list, clean_string




class MissingCrystalError(Exception):
    def __init__(self, id=None):
        if id is None:
            self.message = "Missing Crystal Card information"
        else:
            self.message = f"Missing Crystal Card information for: {id}"

        log("error", self.message)

        super().__init__(self.message)

class SuspiciousCrystalError(MissingCrystalError):
    pass

def parse_crystal_card(file_path) -> dict|None:
    """
    Parse crystal card from .pdb/.cif files.
    :param file_path: File path of file
    :return: Crystal card: { a, b, c,
    alpha, beta, gamma, Z,
    ori1, ori2, ori3,
    scale1, scale2, scale3,
    group_name, group_key }
    """
    log("debug", "Parsing Crystal Card from {}".format(file_path))
    assert os.path.isfile(file_path)
    ext = file_path.split(".")[-1]
    assert "pdb" in ext or "cif" in ext

    if "pdb" in ext:
        cryst1 = None
        origx1 = None
        origx2 = None
        origx3 = None
        scale1 = None
        scale2 = None
        scale3 = None
        with open(file_path) as file:  # Open the file in question
            for line in file:  # Check every line for the keywords
                if "CRYST1" in line: cryst1 = line
                elif "ORIGX1" in line: origx1 = string_to_list(line)
                elif "ORIGX2" in line: origx2 = string_to_list(line)
                elif "ORIGX3" in line: origx3 = string_to_list(line)
                elif "SCALE1" in line: scale1 = string_to_list(line)
                elif "SCALE2" in line: scale2 = string_to_list(line)
                elif "SCALE3" in line: scale3 = string_to_list(line)
        if any([r is None for r in [cryst1, origx1, origx2, origx3, scale1, scale2, scale3]]):
            log("error", "Could not parse Crystal Card from {}".format(file_path))
            return None


        crystal = dict(a=float(clean_string(cryst1[7:16])),
                       b=float(clean_string(cryst1[16:25])),
                       c=float(clean_string(cryst1[25:34])),
                       alpha=float(clean_string(cryst1[34:41])),
                       beta=float(clean_string(cryst1[41:48])),
                       gamma=float(clean_string(cryst1[48:55])),
                       Z=clean_string(cryst1[67:70]),
                       ori1=origx1[1:5], ori2=origx2[1:5], ori3=origx3[1:5],
                       scale1=scale1[1:5], scale2=scale2[1:5], scale3=scale3[1:5])
        crystal["group_name"], crystal["group_key"] = get_space_group(string_to_list(cryst1[55:67])[0:4])

        return crystal

    elif "cif" in ext:
        mmcif = read_mmcif(file_path, subset=["_symmetry", "_cell", "_database_PDB_matrix", "_atom_sites"])
        #print(json.dumps(mmcif.data, indent=4))
        try:
            crystal = dict(

                a=float(mmcif("_cell", "length_a")),
                b=float(mmcif("_cell", "length_b")),
                c=float(mmcif("_cell", "length_c")),
                alpha=float(mmcif("_cell", "angle_alpha")),
                beta=float(mmcif("_cell", "angle_beta")),
                gamma=float(mmcif("_cell", "angle_gamma")),
                Z=float(mmcif("_cell","Z_PDB"))
            )
        except:
            log("error", f"No unit cell found in file: {file_path}")
            raise MissingCrystalError
        try:
            crystal.update(dict(
                group_name=mmcif("_symmetry", "space_group_name_H-M"),
                group_key=int(mmcif("_symmetry", "Int_Tables_number"))
            ))
        except:
            log("warning", f"No symmetry found in file: {file_path}")

        try:
            crystal.update(dict(
                ori1=(
                    float(mmcif("_database_PDB_matrix", "origx[1][1]")),
                    float(mmcif("_database_PDB_matrix", "origx[1][2]")),
                    float(mmcif("_database_PDB_matrix", "origx[1][3]")),
                    float(mmcif("_database_PDB_matrix", "origx_vector[1]"))
                ),
                ori2=(
                    float(mmcif("_database_PDB_matrix", "origx[2][1]")),
                    float(mmcif("_database_PDB_matrix", "origx[2][2]")),
                    float(mmcif("_database_PDB_matrix", "origx[2][3]")),
                    float(mmcif("_database_PDB_matrix", "origx_vector[2]"))
                ),
                ori3=(
                    float(mmcif("_database_PDB_matrix", "origx[3][1]")),
                    float(mmcif("_database_PDB_matrix", "origx[3][2]")),
                    float(mmcif("_database_PDB_matrix", "origx[3][3]")),
                    float(mmcif("_database_PDB_matrix", "origx_vector[3]"))
                )
            ))
        except TypeError:
            log("warning", "Error parsing crystal origin")
        try:
            crystal.update(dict(
                scale1=(
                    float(mmcif("_atom_sites", "fract_transf_matrix[1][1]")),
                    float(mmcif("_atom_sites", "fract_transf_matrix[1][2]")),
                    float(mmcif("_atom_sites", "fract_transf_matrix[1][3]")),
                    float(mmcif("_atom_sites", "fract_transf_vector[1]"))
                ),
                scale2=(
                    float(mmcif("_atom_sites", "fract_transf_matrix[2][1]")),
                    float(mmcif("_atom_sites", "fract_transf_matrix[2][2]")),
                    float(mmcif("_atom_sites", "fract_transf_matrix[2][3]")),
                    float(mmcif("_atom_sites", "fract_transf_vector[2]"))
                ),
                scale3=(
                    float(mmcif("_atom_sites", "fract_transf_matrix[3][1]")),
                    float(mmcif("_atom_sites", "fract_transf_matrix[3][2]")),
                    float(mmcif("_atom_sites", "fract_transf_matrix[3][3]")),
                    float(mmcif("_atom_sites", "fract_transf_vector[3]"))
                )
            ))
        except TypeError:
            log("warning", "Error parsing crystal scale")

        if crystal["a"] == crystal["b"] == crystal["c"]:
            if crystal["alpha"] == crystal["beta"] == crystal["gamma"]:
                if crystal["Z"] == 1.:
                    if "P" in crystal["group_name"] and "1" in crystal["group_name"]:
                        if crystal["ori1"][0] + crystal["ori2"][1] + crystal["ori3"][2] == 3.:
                            if crystal["scale1"][0] + crystal["scale2"][1] + crystal["scale3"][2] == 3.:
                                raise SuspiciousCrystalError
        #print(json.dumps(crystal, indent=4))
        return crystal


    log("error", "Could not parse Crystal Card from {}".format(file_path))
    return None


def get_space_group(raw_group:str) -> list[str|int]:
    """
    Parse space group from PDB file CRYST1 line.
    :param raw_group: Raw space group name
    :return: [space group name, space group key]
    """
    new_group = None
    new_key = None
    raw_group = " ".join(raw_group).strip()
    log("debug", raw_group, end=" -> ")
    for key, group in dictio_space_groups.items():
        if raw_group == group["symbol"]:
            new_group = str(group["symbol"])
            new_key = key
            break
    if new_group is None or new_key is None:
        log("error", "Missing crystal group")
        return None
    else:
        log("debug", new_group)

    return [new_group, new_key]


def get_cell_dim(card:dict) -> list:
    """
    Parse unit cell parameters from crystal card.
    :param card: Crystal card
    :return: [A, B, C, alpha, beta, gamma]
    """
    return [card["a"], card["b"], card["c"], card["alpha"], card["beta"], card["gamma"]]



def calculate_parameters(card:dict) -> dict:
    """
    Calculate parameters from crystal card.
    :param card: Crystal card
    :return: Parameters dictionary
    """
    if card is None:
        log("warning", "No card loaded")
        return None
    cell_dim = get_cell_dim(card)
    parameters = {}
    parameters["A"] = A = float(cell_dim[0])
    parameters["B"] = B = float(cell_dim[1])
    parameters["C"] = C = float(cell_dim[2])
    parameters["alphaDeg"] = alphaDeg = float(cell_dim[3])
    parameters["betaDeg"] = betaDeg = float(cell_dim[4])
    parameters["gammaDeg"] = gammaDeg = float(cell_dim[5])
    parameters["alpha"] = alpha = (alphaDeg * 2 * np.pi) / 360
    parameters["beta"] = beta = (betaDeg * 2 * np.pi) / 360
    parameters["gamma"] = gamma = (gammaDeg * 2 * np.pi) / 360
    parameters["c_a"] = c_a = np.cos(alpha)
    parameters["c_b"] = c_b = np.cos(beta)
    parameters["c_g"] = c_g = np.cos(gamma)
    parameters["s_g"] = s_g = np.sin(gamma)
    parameters["q"] = q = np.sqrt(1 + 2 * c_a * c_b * c_g - c_a ** 2 - c_b ** 2 - c_g ** 2)
    parameters["uu"] = uu = s_g / (q * C)
    parameters["vv"] = vv = (c_b * c_g - c_a) / (q * B * s_g)
    parameters["uuy"] = uuy = 1 / (B * s_g)
    parameters["vvz"] = vvz = -1 * (c_g / (A * s_g))
    parameters["uuz"] = uuz = (c_a * c_g - c_b) / (q * A * s_g)
    parameters["vvy"] = vvy = 1 / A

    return parameters


