import os, sys, math, json

import numpy as np
from ..utilities.kdtree import KDT

from ..base import PseudoAtom
from ..utilities import log

atomic_radii = {
    "default": {
        "H": 1.200,
        "HE": 1.400,
        "C": 1.700,
        "N": 1.550,
        "O": 1.520,
        "F": 1.470,
        "NA": 2.270,
        "MG": 1.730,
        "P": 1.800,
        "S": 1.800,
        "CL": 1.750,
        "K": 2.750,
        "CA": 2.310,
        "NI": 1.630,
        "CU": 1.400,
        "ZN": 1.390,
        "SE": 1.900,
        "BR": 1.850,
        "CD": 1.580,
        "I": 1.980,
        "HG": 1.550,
        "other": 2.000,
    }
}

residue_sasas = {
    "ALA": 129.0,
    "ARG": 274.0,
    "ASN": 195.0,
    "ASP": 193.0,
    "CYS": 167.0,
    "GLU": 223.0,
    "GLN": 225.0,
    "GLY": 104.0,
    "HIS": 224.0,
    "ILE": 197.0,
    "LEU": 201.0,
    "LYS": 236.0,
    "MET": 224.0,
    "PHE": 240.0,
    "PRO": 159.0,
    "SER": 155.0,
    "THR": 172.0,
    "TRP": 285.0,
    "TYR": 263.0,
    "VAL": 174.0,
}



class SASA(object):
    def __init__(self, ball_radius=1.40, n_points=100, radii_dict="default", **kwargs):

        assert ball_radius > 0
        assert n_points > 1

        self.ball_radius = ball_radius
        self.n_points = n_points
        global atomic_radii
        self.radii_dict_name = radii_dict
        self.radii_dict = atomic_radii[radii_dict]

        self._sphere = self._compute_sphere()

    def _compute_sphere(self):
        n = self.n_points

        dl = np.pi * (3 - 5 ** 0.5)
        dz = 2.0 / n

        longitude = 0
        z = 1 - dz / 2

        coords = np.zeros((n, 3), dtype=np.float32)
        for k in range(n):
            r = (1 - z * z) ** 0.5
            coords[k, 0] = math.cos(longitude) * r
            coords[k, 1] = math.sin(longitude) * r
            coords[k, 2] = z
            z -= dz
            longitude += dl

        return coords

    def compute(self, entity, targets=None, save_sasas=None, force=False, quiet=False, **kwargs):
        if not quiet:
            log(2, "Computing SASA...")

        if getattr(entity, "_kdtrees", {}).get("atom", None) is not None and not force:
            if not quiet:
                log(3, "Recovering saved KDTree (atom)")
            kdt = entity._kdtrees["atom"]
        else:
            kdt = KDT(entity, quiet=quiet, kdtree_name="atom", **kwargs)

        radii_list = []
        n_atoms = 0
        for a in kdt.atoms:
            n_atoms += 1
            if a.element in self.radii_dict:
                radii_list.append(self.radii_dict[a.element])
            else:
                radii_list.append(self.radii_dict["other"])

        radii = np.array(radii_list, dtype=np.float64)


        radii += self.ball_radius
        twice_maxradii = np.max(radii) * 2

        # ptset = set(range(self.n_points))

        # print(ptset)

        if targets is None:
            if save_sasas is None:
                save_sasas = True
            targets = kdt.atoms
        else:
            targets = [PseudoAtom(c) if not isinstance(c, PseudoAtom) else c for c in targets]
            if save_sasas is None:
                save_sasas = False

        asa_array = np.zeros((len(targets), 1), dtype=np.int64)
        target_radii = []
        for a in targets:
            if a.element in self.radii_dict:
                target_radii.append(self.radii_dict[a.element])
            else:
                target_radii.append(self.radii_dict["other"])
        target_radii = np.array(target_radii, dtype=np.float64)


        for i, target in enumerate(targets):
            log(3, f"{i}/{len(targets)}", end="\r")
            # exposed_points = ptset.copy()

            if isinstance(target, PseudoAtom):
                target = target.coord

            i_radii = radii[i]
            s_on_i = (np.array(self._sphere, copy=True) * i_radii) + np.array(target)

            sphere_kdt = KDT(s_on_i, leaf_size=10, quiet=True)

            i_neighbours = kdt.of(coords=s_on_i, radius=twice_maxradii, distances=False, unique=True)

            # print(i, len(i_neighbours))
            neighbour_radii = np.array([radii_list[j] for j in i_neighbours])
            neighbour_coords = np.array([kdt.coords[j] for j in i_neighbours])

            overlap_indexes = sphere_kdt.of(neighbour_coords, radius=neighbour_radii, unique=True)

            # print(i, len(overlap_indexes))

            asa_array[i] = self.n_points - len(overlap_indexes)
            #print(asa_array[i], self.n_points, len(overlap_indexes), self.n_points - len(overlap_indexes))


        f = target_radii * target_radii * (4 * np.pi / self.n_points)
        asa_array = asa_array * f[:, np.newaxis]

        if save_sasas:
            for atom, asa in zip(targets, asa_array):
                if isinstance(atom, PseudoAtom):
                    atom.set_misc("SASA", float(asa[0]))
            entity.set_flag("sasa_calculated", True)
            entity.data["SASA"]["ball_radius"] = self.ball_radius
            entity.data["SASA"]["n_points"] = self.n_points
            entity.data["SASA"]["radii_dict"] = self.radii_dict_name
            entity.data["SASA"]["average"] = None
        return asa_array








