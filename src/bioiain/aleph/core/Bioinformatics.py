#! /usr/bin/env python
# -*- coding: utf-8 -*-

#future imports
from __future__ import print_function
from __future__ import division

#from __future__ import unicode_literals

#System imports
import sys
import os
import copy
import re
import shutil
import datetime
from collections import defaultdict
from collections import OrderedDict
import heapq

# Python standard modules imports
import pickle
import operator
import itertools
import traceback
import subprocess

import io
import gzip
import shutil

import warnings
warnings.filterwarnings("ignore")

#Scientific and numerical imports
import Bio.PDB
from Bio.PDB import Entity
from Bio import SeqIO
import numpy
import scipy
import scipy.stats
#import csb.bio.utils
import scipy.cluster.vq

#Other imports
from termcolor import colored
import igraph
#import igraph.vendor.texttable

#ARCIMBOLDO_FULL imports
from . import ALEPH
if os.path.basename(ALEPH.__file__)=="__init__.py":
    import ALEPH.aleph.core.ALEPH as ALEPH

try:
    import urllib
    from builtins import range
    from builtins import str
except ImportError:
    import urllib2

from . import SystemUtility

#TODO: This function is needed by getSuperimp but must be deleted with it
SUFRAGLENGTH = 3

AADICMAP = OrderedDict([('ALA', 'A'), ('CYS', 'C'), ('CSO', 'C'), ('OCS', 'C'), ('ASP', 'D'), ('GLU', 'E'), ('PHE', 'F'), ('GLY', 'G'),
                         ('HIS', 'H'), ('ILE', 'I'), ('LYS', 'K'), ('LEU', 'L'), ('MET', 'M'), ('MSE', 'M'), ('SAM', 'M'), ('ASN', 'N'),
                         ('PRO', 'P'), ('GLN', 'Q'), ('ARG', 'R'), ('SER', 'S'), ('THR', 'T'), ('KCX', 'K'), ('VAL', 'V'), ('TRP', 'W'),
                         ('TYR', 'Y'), ('UNK', '-'), ('SEC', 'U'), ('PYL', 'O'), ('ASX', 'B'), ('GLX', 'Z'), ('XLE', 'J'), ('FOL', 'X'),
                         ('DAL', 'A'), ('XAA', 'X'), ('ME0', 'M'), ('CGU', '-')])

ATOAAAMAP =  dict((v,k) for k,v in list(AADICMAP.items())[::-1])


AADICMAP = defaultdict(lambda: "X", **AADICMAP)

AAList = AADICMAP.keys()

AALISTOL = AADICMAP.values()

list_id = 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890'  # All possible chain ids for a PDB

def parse_property_table(table_string):
    value_dict = {}
    for line in table_string.splitlines():
        line = line.strip()
        if not line:
            continue
        fields = line.split(" ")
        fields = [f for f in fields if len(f.strip()) > 0]
        assert len(fields) >= 2
        value, letter = fields[:2]
        assert letter not in value_dict, "Repeated amino acid " + line
        value_dict[letter] = float(value)
    return value_dict

"""
Amino acids property tables copied from CRASP website

###
# Values copied from:
# "Solvent accessibility of AA in known protein structures"
# http://prowl.rockefeller.edu/aainfo/access.htm
###
"""

hydropathy = parse_property_table("""
1.80000 A ALA
-4.5000 R ARG
-3.5000 N ASN
-3.5000 D ASP
2.50000 C CYS
-3.5000 Q GLN
-3.5000 E GLU
-0.4000 G GLY
-3.2000 H HIS
4.50000 I ILE
3.80000 L LEU
-3.9000 K LYS
1.90000 M MET
2.80000 F PHE
-1.6000 P PRO
-0.8000 S SER
-0.7000 T THR
-0.9000 W TRP
-1.3000 Y TYR
4.20000 V VAL
""")
volume = parse_property_table("""
91.5000 A ALA
202.0000 R ARG
135.2000 N ASN
124.5000 D ASP
118.0000 C CYS
161.1000 Q GLN
155.1000 E GLU
66.40000 G GLY
167.3000 H HIS
168.8000 I ILE
167.9000 L LEU
171.3000 K LYS
170.8000 M MET
203.4000 F PHE
129.3000 P PRO
99.10000 S SER
122.1000 T THR
237.6000 W TRP
203.6000 Y TYR
141.7000 V VAL
""")
polarity = parse_property_table("""
0.0000 A ALA
52.000 R ARG
3.3800 N ASN
40.700 D ASP
1.4800 C CYS
3.5300 Q GLN
49.910 E GLU
0.0000 G GLY
51.600 H HIS
0.1500 I ILE
0.4500 L LEU
49.500 K LYS
1.4300 M MET
0.3500 F PHE
1.5800 P PRO
1.6700 S SER
1.6600 T THR
2.1000 W TRP
1.6100 Y TYR
0.1300 V VAL
""")
pK_side_chain = parse_property_table("""
0.0000 A ALA
12.480 R ARG
0.0000 N ASN
3.6500 D ASP
8.1800 C CYS
0.0000 Q GLN
4.2500 E GLU
0.0000 G GLY
6.0000 H HIS
0.0000 I ILE
0.0000 L LEU
10.530 K LYS
0.0000 M MET
0.0000 F PHE
0.0000 P PRO
0.0000 S SER
0.0000 T THR
0.0000 W TRP
10.700 Y TYR
0.0000 V VAL
""")
prct_exposed_residues = parse_property_table("""
15.0000 A ALA
67.0000 R ARG
49.0000 N ASN
50.0000 D ASP
5.00000 C CYS
56.0000 Q GLN
55.0000 E GLU
10.0000 G GLY
34.0000 H HIS
13.0000 I ILE
16.0000 L LEU
85.0000 K LYS
20.0000 M MET
10.0000 F PHE
45.0000 P PRO
32.0000 S SER
32.0000 T THR
17.0000 W TRP
41.0000 Y TYR
14.0000 V VAL
""")
hydrophilicity = parse_property_table("""
-0.5000 A ALA
3.00000 R ARG
0.20000 N ASN
3.00000 D ASP
-1.0000 C CYS
0.20000 Q GLN
3.00000 E GLU
0.00000 G GLY
-0.5000 H HIS
-1.8000 I ILE
-1.8000 L LEU
3.00000 K LYS
-1.3000 M MET
-2.5000 F PHE
0.00000 P PRO
0.30000 S SER
-0.4000 T THR
-3.4000 W TRP
-2.3000 Y TYR
-1.5000 V VAL
""")
accessible_surface_area = parse_property_table("""
27.8000 A ALA
94.7000 R ARG
60.1000 N ASN
60.6000 D ASP
15.5000 C CYS
68.7000 Q GLN
68.2000 E GLU
24.5000 G GLY
50.7000 H HIS
22.8000 I ILE
27.6000 L LEU
103.000 K LYS
33.5000 M MET
25.5000 F PHE
51.5000 P PRO
42.0000 S SER
45.0000 T THR
34.7000 W TRP
55.2000 Y TYR
23.7000 V VAL
""")
local_flexibility = parse_property_table("""
705.42000 A ALA
1484.2800 R ARG
513.46010 N ASN
34.960000 D ASP
2412.5601 C CYS
1087.8300 Q GLN
1158.6600 E GLU
33.180000 G GLY
1637.1300 H HIS
5979.3701 I ILE
4985.7300 L LEU
699.69000 K LYS
4491.6602 M MET
5203.8599 F PHE
431.96000 P PRO
174.76000 S SER
601.88000 T THR
6374.0698 W TRP
4291.1001 Y TYR
4474.4199 V VAL
""")
accessible_surface_area_folded = parse_property_table("""
31.5000 A ALA
93.8000 R ARG
62.2000 N ASN
60.9000 D ASP
13.9000 C CYS
74.0000 Q GLN
72.3000 E GLU
25.2000 G GLY
46.7000 H HIS
23.0000 I ILE
29.0000 L LEU
110.300 K LYS
30.5000 M MET
28.7000 F PHE
53.7000 P PRO
44.2000 S SER
46.0000 T THR
41.7000 W TRP
59.1000 Y TYR
23.5000 V VAL
""")
refractivity = parse_property_table("""
4.34000 A ALA
26.6600 R ARG
13.2800 N ASN
12.0000 D ASP
35.7700 C CYS
17.5600 Q GLN
17.2600 E GLU
0.00000 G GLY
21.8100 H HIS
19.0600 I ILE
18.7800 L LEU
21.2900 K LYS
21.6400 M MET
29.4000 F PHE
10.9300 P PRO
6.35000 S SER
11.0100 T THR
42.5300 W TRP
31.5300 Y TYR
13.9200 V VAL
""")
mass = parse_property_table("""
70.079 A ALA
156.188 R ARG
114.104 N ASN
115.089 D ASP
103.144 C CYS
128.131 Q GLN
129.116 E GLU
57.052 G GLY
137.142 H HIS
113.160 I ILE
113.160 L LEU
128.174 K LYS
131.198 M MET
147.177 F PHE
97.177 P PRO
87.078 S SER
101.105 T THR
186.213 W TRP
163.170 Y TYR
99.133 V VAL
""")
solvent_exposed_area = dict(
    S=0.70,
    T=0.71,
    A=0.48,
    G=0.51,
    P=0.78,
    C=0.32,
    D=0.81,
    E=0.93,
    Q=0.81,
    N=0.82,
    L=0.41,
    I=0.39,
    V=0.40,
    M=0.44,
    F=0.42,
    Y=0.67,
    W=0.49,
    K=0.93,
    R=0.84,
    H=0.66,
)

atom_weights = {
    'H': 1.00794,
    'He': 4.002602,
    'Li': 6.941,
    'Be': 9.012182,
    'B': 10.811,
    'C': 12.0107,
    'N': 14.0067,
    'O': 15.9994,
    'F': 18.9984032,
    'Ne': 20.1797,
    'Na': 22.989770,
    'Mg': 24.3050,
    'Al': 26.981538,
    'Si': 28.0855,
    'P': 30.973761,
    'S': 32.065,
    'Cl': 35.453,
    'Ar': 39.948,
    'K': 39.0983,
    'Ca': 40.078,
    'Sc': 44.955910,
    'Ti': 47.867,
    'V': 50.9415,
    'Cr': 51.9961,
    'Mn': 54.938049,
    'Fe': 55.845,
    'Co': 58.933200,
    'Ni': 58.6934,
    'Cu': 63.546,
    'Zn': 65.39,
    'Ga': 69.723,
    'Ge': 72.64,
    'As': 74.92160,
    'Se': 78.96,
    'Br': 79.904,
    'Kr': 83.80,
    'Rb': 85.4678,
    'Sr': 87.62,
    'Y': 88.90585,
    'Zr': 91.224,
    'Nb': 92.90638,
    'Mo': 95.94,
    'Tc': 98.0,
    'Ru': 101.07,
    'Rh': 102.90550,
    'Pd': 106.42,
    'Ag': 107.8682,
    'Cd': 112.411,
    'In': 114.818,
    'Sn': 118.710,
    'Sb': 121.760,
    'Te': 127.60,
    'I': 126.90447,
    'Xe': 131.293,
    'Cs': 132.90545,
    'Ba': 137.327,
    'La': 138.9055,
    'Ce': 140.116,
    'Pr': 140.90765,
    'Nd': 144.24,
    'Pm': 145.0,
    'Sm': 150.36,
    'Eu': 151.964,
    'Gd': 157.25,
    'Tb': 158.92534,
    'Dy': 162.50,
    'Ho': 164.93032,
    'Er': 167.259,
    'Tm': 168.93421,
    'Yb': 173.04,
    'Lu': 174.967,
    'Hf': 178.49,
    'Ta': 180.9479,
    'W': 183.84,
    'Re': 186.207,
    'Os': 190.23,
    'Ir': 192.217,
    'Pt': 195.078,
    'Au': 196.96655,
    'Hg': 200.59,
    'Tl': 204.3833,
    'Pb': 207.2,
    'Bi': 208.98038,
    'Po': 208.98,
    'At': 209.99,
    'Rn': 222.02,
    'Fr': 223.02,
    'Ra': 226.03,
    'Ac': 227.03,
    'Th': 232.0381,
    'Pa': 231.03588,
    'U': 238.02891,
    'Np': 237.05,
    'Pu': 244.06,
    'Am': 243.06,
    'Cm': 247.07,
    'Bk': 247.07,
    'Cf': 251.08,
    'Es': 252.08,
    'Fm': 257.10,
    'Md': 258.10,
    'No': 259.10,
    'Lr': 262.11,
    'Rf': 261.11,
    'Db': 262.11,
    'Sg': 266.12,
    'Bh': 264.12,
    'Hs': 269.13,
    'Mt': 268.14,
}


def center_of_mass(entity, geometric=False):
    """
    Returns gravitic [default] or geometric center of mass of an Entity.
    Geometric assumes all masses are equal (geometric=True)
    """
    # Structure, Model, Chain, Residue
    if isinstance(entity, Entity.Entity):
        atom_list = entity.get_atoms()
    # List of Atoms
    elif hasattr(entity, '__iter__') and [x for x in entity if x.level == 'A']:
        atom_list = entity
    else:  # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                         "Structure, Model, Chain, Residue, list of Atoms.")

    class COM:
        def __init__(self, coord):
            self.coord = coord

    positions = [[], [], []]  # [ [X1, X2, ..] , [Y1, Y2, ...] , [Z1, Z2, ...] ]
    masses = []

    for atom in atom_list:
        try:
            atom.mass = atom_weights[atom.element.capitalize()]
        except:
            atom.mass = 1.0

        masses.append(atom.mass)

        for i, coord in enumerate(atom.coord.tolist()):
            positions[i].append(coord)

    # If there is a single atom with undefined mass complain loudly.
    if 'ukn' in set(masses) and not geometric:
        raise ValueError("Some Atoms don't have an element assigned.\n"
                         "Try adding them manually or calculate the geometrical center of mass instead.")

    if geometric:
        com = COM([sum(coord_list) / len(masses) for coord_list in positions])
        return com
    else:
        w_pos = [[], [], []]
        for atom_index, atom_mass in enumerate(masses):
            w_pos[0].append(positions[0][atom_index] * atom_mass)
            w_pos[1].append(positions[1][atom_index] * atom_mass)
            w_pos[2].append(positions[2][atom_index] * atom_mass)
        com = COM([sum(coord_list) / sum(masses) for coord_list in w_pos])
        return com


def calculate_shape_param(structure):
    """
    Calculates the gyration tensor of a structure.
    Returns a tuple containing shape parameters:

      ((a,b,c), rg, A, S)

      (a,b,c) - dimensions of the semi-axis of the ellipsoid
      rg    - radius of gyration of the structure
      A     - sphericity value
      S     - anisotropy value

      References:
           1.  Rawat N, Biswas P (2011) Shape, flexibility and packing of proteins and nucleic acids in complexes. Phys Chem Chem Phys 13:9632-9643
           2.  Thirumalai D (2004) Asymmetry in the Shapes of Folded and Denatured States of Proteinss - The Journal of Physical Chemistry B
           3.  Vondrasek J (2011) Gyration- and Inertia-Tensor-Based Collective Coordinates for Metadynamics. Application on the Conformational Behavior of Polyalanine Peptides and Trp-Cage Folding - The Journal of Physical Chemistry A
    """

    com = center_of_mass(structure, True)
    cx, cy, cz = com.coord

    n_atoms = 0
    tensor_xx, tensor_xy, tensor_xz = 0, 0, 0
    tensor_yx, tensor_yy, tensor_yz = 0, 0, 0
    tensor_zx, tensor_zy, tensor_zz = 0, 0, 0

    if isinstance(structure, Entity.Entity):
        atom_list = structure.get_atoms()
    # List of Atoms
    elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
        atom_list = structure
    else:  # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                         "Structure, Model, Chain, Residue, list of Atoms.")

    for atom in atom_list:
        ax, ay, az = atom.coord
        tensor_xx += (ax - cx) * (ax - cx)
        tensor_yx += (ax - cx) * (ay - cy)
        tensor_xz += (ax - cx) * (az - cz)
        tensor_yy += (ay - cy) * (ay - cy)
        tensor_yz += (ay - cy) * (az - cz)
        tensor_zz += (az - cz) * (az - cz)
        n_atoms += 1

    gy_tensor = numpy.mat(
        [[tensor_xx, tensor_yx, tensor_xz], [tensor_yx, tensor_yy, tensor_yz], [tensor_xz, tensor_yz, tensor_zz]])
    gy_tensor = (1.0 / n_atoms) * gy_tensor

    D, V = numpy.linalg.eig(gy_tensor)
    [a, b, c] = sorted(numpy.sqrt(D))  # lengths of the ellipsoid semi-axis
    rg = numpy.sqrt(sum(D))

    l = numpy.average([D[0], D[1], D[2]])
    A = (((D[0] - l) ** 2 + (D[1] - l) ** 2 + (D[2] - l) ** 2) / l ** 2) * 1 / 6
    S = (((D[0] - l) * (D[1] - l) * (D[2] - l)) / l ** 3)

    return ((a * 2, b * 2, c * 2), rg, A, S)

    # print "%s" % '#Dimensions(a,b,c) #Rg #Anisotropy'
    # print "%.2f" % round(a,2), round(b,2), round(c,2) , round(rg,2) , round(A,2)


def calculate_moment_of_intertia_tensor(structure):
    """
    Calculates the moment of inertia tensor from the molecule.
    Returns a numpy matrix.
    """

    if isinstance(structure, Entity.Entity):
        atom_list = structure.get_atoms()
        # List of Atoms
    elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
        atom_list = structure
    else:  # Some other weirdo object
        raise ValueError("Center of Mass can only be calculated from the following objects:\n"
                         "Structure, Model, Chain, Residue, list of Atoms.")
    s_mass = 0.0
    for atom in atom_list:
        atom.mass = atom_weights[atom.element.capitalize()]
        s_mass += atom.mass

    com = center_of_mass(structure, False)
    cx, cy, cz = com.coord

    n_atoms = 0
    tensor_xx, tensor_xy, tensor_xz = 0, 0, 0
    tensor_yx, tensor_yy, tensor_yz = 0, 0, 0
    tensor_zx, tensor_zy, tensor_zz = 0, 0, 0
    # s_mass = sum([a.mass for a in atom_list])

    if isinstance(structure, Entity.Entity):
        atom_list = structure.get_atoms()
    elif hasattr(structure, '__iter__') and [x for x in structure if x.level == 'A']:
        atom_list = structure

    for atom in atom_list:
        ax, ay, az = atom.coord
        tensor_xx += ((ay - cy) ** 2 + (az - cz) ** 2) * atom_weights[atom.element.capitalize()]
        tensor_xy += -1 * (ax - cx) * (ay - cy) * atom_weights[atom.element.capitalize()]
        tensor_xz += -1 * (ax - cx) * (az - cz) * atom_weights[atom.element.capitalize()]
        tensor_yy += ((ax - cx) ** 2 + (az - cz) ** 2) * atom_weights[atom.element.capitalize()]
        tensor_yz += -1 * (ay - cy) * (az - cz) * atom_weights[atom.element.capitalize()]
        tensor_zz += ((ax - cx) ** 2 + (ay - cy) ** 2) * atom_weights[atom.element.capitalize()]

    in_tensor = numpy.mat([[tensor_xx, tensor_xy, tensor_xz], [tensor_xy, tensor_yy, tensor_yz], [tensor_xz,
                                                                                            tensor_yz, tensor_zz]])
    D, V = numpy.linalg.eig(in_tensor)

    a = numpy.sqrt((5 / (2 * s_mass)) * (D[0] - D[1] + D[2]))
    b = numpy.sqrt((5 / (2 * s_mass)) * (D[2] - D[0] + D[1]))
    c = numpy.sqrt((5 / (2 * s_mass)) * (D[1] - D[2] + D[0]))
    return sorted([a, b, c]), D, V


def set_occupancy_to_zero_for_outliers(structure,model,dizio_resi):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure: Bio.PDB.Structure
    :param model: int
    :param dizio_resi: dict
    :return: Bio.PDB.Structure
    """
    for resi in get_list_of_residues(structure,model):
        if resi.get_full_id()[1:4] not in [ren[1][1:4] for ren in dizio_resi.values()]:
            for atom in resi:
                atom.set_occupancy(0.0)
    return structure

def trim_ensemble_to_distance(pdb,codec,distance):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param pdb:
    :param codec:
    :param distance:
    :return:
    """
    if (isinstance(pdb,str) and os.path.exists(pdb)) or isinstance(pdb,io.StringIO):
        if isinstance(pdb,io.StringIO):
            pdb.seek(0)
        structure1 = get_structure("a",pdb)
        list_CA1 = [residue['CA'] for residue in Bio.PDB.Selection.unfold_entities(structure1,'R') if residue.has_id("CA")]
    else:
        list_CA1 = pdb

    distance_hash = {ca1.get_full_id(): [(ca1-ca2,ca2.get_full_id(),ca1,ca2) for ca2 in list_CA1 if ca1.get_full_id()[1]!=ca2.get_full_id()[1]] for ca1 in list_CA1}
    eliminate = []
    pairs = {}
    for key in distance_hash:
        lista = sorted(distance_hash[key], key=lambda x: (x[0],x[1][1]))
        models = []
        zer = []
        for one in lista:
            if one[1][1] not in models and one[0]<=distance and one[0]>0:
                zer.append(one)
                models.append((one[1][1],one[1]))

        if key[1] not in pairs:
            pairs[key[1]] = {m[0]:[m[1]] for m in models}
        else:
            for tk in models:
                if tk[0] not in pairs[key[1]]:
                    pairs[key[1]][tk[0]] = [tk[1]]
                else:
                    pairs[key[1]][tk[0]].append(tk[1])

        if len(zer) > 0:
            distance_hash[key] = zer
        else:
            eliminate.append(key)

    for key in eliminate:
        del distance_hash[key]

    atoms = set([])
    for cosmo in distance_hash.keys():
        atoms.add((cosmo, distance_hash[cosmo][0][2]))
        for idr in distance_hash[cosmo]:
            atoms.add((idr[1], idr[3]))
    atoms = sorted(list(atoms))
    models = {}
    for atm in atoms:
        if atm[0][1] not in models:
            models[atm[0][1]] = [atm[1]]
        else:
            models[atm[0][1]].append(atm[1])

    size_frag = 2
    for model in models:
        models[model] = sorted(list(set(models[model])),key=lambda x: x.get_full_id())
        remove = []
        #print("This is model",model)
        cont = 0
        frags = set([])
        for e in range(len(models[model])-1):
            frags.add(e)
            if check_continuity(models[model][e].get_parent(), models[model][e + 1].get_parent(), verbose=False):
                cont += 1
                frags.add(e+1)
            else:
                if cont < size_frag:
                    for c in frags:
                        #print("Removing..",models[model][c].get_full_id())
                        remove.append(c)
                frags = set([])
                cont = 0

        remove = sorted(remove,reverse=True)
        for rem in remove:
            del models[model][rem]


    #quest = []
    kk_to_del = []
    for key,value in pairs.items():
        #print("----",key,value)
        key_to_del = []
        for v in value.keys():
            value[v] = sum([1 for t in value[v] if t in [s.get_full_id() for s in models[v]]])
            if value[v] == 0:
                key_to_del.append(v)

        for v in key_to_del:
            del value[v]

        #print("++++",key,value)
        pairs[key] = value
        if len(value.keys()) == 0:
            kk_to_del.append(key)

        #TODO: here I can force to check a minimum amount of number of contacts to allow insertion in the quest list
        #quest.append([key]+list(value.keys()))

    for v in kk_to_del:
        del pairs[v]

    all_reductions = []
    for p in sorted(pairs.keys(), key=lambda x: len(pairs[x].keys()), reverse=True):
        queue = [(p,0)]
        for t in sorted(pairs[p].keys(), key=lambda x: pairs[p][x], reverse=True):
            insert = True
            for c in queue:
                if t not in pairs[c[0]]:
                    insert = False
                    break
            if insert:
                queue.append((t, pairs[p][t]))
        print("Reduction:",queue)
        all_reductions.append(queue)

    if len(all_reductions) <= 0:
        return ""

    result = set(sorted(all_reductions, key=lambda x: (len(x), sum([o[1] for o in x])), reverse=True)[0])

    # quest = sorted(quest,key=lambda x: len(x), reverse=True)
    #
    # result = set(quest[0])
    # visited = [quest[0][0]]
    # for s in quest[1:]:
    #     tr = result.copy()
    #     result.intersection_update(s)
    #     visited.append(s[0])
    #     print("Exploring intersection:",result)
    #     if all([pr in visited for pr in result]):
    #         break
    #     if len(result) < 4:
    #         result = tr
    #         print("going back to",result)
    #         break
    print("Final models selected",result)
    result = set([o[0] for o in result])

    #distance_hash2 = {key: sorted(distance_hash[key], key=lambda x: (x[0],x[1][1]))[0] for key in distance_hash if (sorted(distance_hash[key], key=lambda x: (x[0],x[1][1]))[0][0] <= distance)}

    pdball = ""
    for model in sorted(models.keys()):
        if model in result:
            pdball += "MODEL " + str(model) + "\n"
            pdball += "REMARK TITLE " + codec[model] + "\n"
            pdball += get_pdb_from_list_of_atoms(models[model], renumber=False, uniqueChain=False, chainId="A", chainFragment=False, diffchain=None, polyala=False, maintainCys=False, normalize=False, sort_reference=True)[0]
            pdball += "ENDMDL" + "\n"

    return pdball

def distance_between_matrices(alist,blist,are_cosins=False,verbose=False, weights_list=None, stat_test=True):
    #numpy.seterr(all='raise')
    try:
        a = numpy.concatenate([cc.flatten() for cc in alist])
        b = numpy.concatenate([dd.flatten() for dd in blist])

        if weights_list is not None:
            weights = numpy.concatenate([ee.flatten() for ee in weights_list])
        else:
            weights = None
        e = None
        if stat_test and a.shape[0]>20:
        #    #if not t_student(a,b,weights=weights, p_value=0.05): return 0.0 #Not comparable too different
            e = mann_whitney(a,b, are_cosins=are_cosins, weights=None, p_value=0.1, verbose=False)

        sigma = 1
        if a.shape[0]<20:  sigma = 2

        if not are_cosins:
            ddd = a - b
            max_one = a.copy()
            max_one[numpy.where(ddd < 0)] = b[numpy.where(ddd < 0)]
            Z = numpy.abs(ddd) / max_one

            if weights is not None and weights.shape == Z.shape:
                mean = numpy.average(Z,weights=weights)
                std = numpy.sqrt(numpy.average((Z-mean)**2, weights=weights))
            else:
                mean = numpy.mean(Z)
                std = numpy.std(Z)
            diff = 1.0 - (mean+(sigma*std))
            if e is not None :
                if e: diff += 0.1
                elif diff<0.8: diff -= 0.2 #0.15
            #if a.shape[0]<20:  diff -= 0.1
        else:
            #NOTE: The current defaults values are already normalized in (0, 1) from the uniform distributed degree angles

            #NOTE: The following normalization is to trasform cosins (-1, 1) to (0,1). Trouble is that cosin is not uniformely distributed
            # a = (((a+1.0) * 2) / 2.0) / 2.0
            # b = (((b+1.0) * 2) / 2.0) / 2.0
            #NOTE: The following normalization is to transform angles in radians in (0, 1). Tried but It does not work as good as with degrees
            #a /= numpy.pi
            #b /= numpy.pi

            ddd = a - b
            Z = numpy.abs(ddd) #/ max_one #/ 2.0 #/ (numpy.pi) #numpy.pi or 2*numpy.pi???
            if weights is not None and weights.shape == Z.shape:
                mean = numpy.average(Z, weights=weights)
                std = numpy.sqrt(numpy.average((Z - mean) ** 2, weights=weights))
            else:
                mean = numpy.mean(Z)
                std = numpy.std(Z)
            diff = 1.0 - (mean+(sigma*std))
            if e is not None:
                if e: diff += 0.1
                elif diff<0.8: diff -= 0.2 #0.15
            #if a.shape[0]<20:  diff -= 0.1

        #verbose=True
        # if weights is not None:
        #     print(list(Z))
        #     print(list(weights))
        if verbose: #and diff>0.3 and not are_cosins: #and diff>0.10:
            print("A", a)
            print("B", b)
            print("a-b", abs(a - b))
            print("Z",Z)
            #print("maxone", max_one)
            print("ANGLE:" if are_cosins else "DIST:","Mean",numpy.mean(Z),"STD",numpy.std(Z),"MEAN + STD",(numpy.mean(Z)+numpy.std(Z)),
                  "MeanW",mean,"STDW",std,"MEANW + STDW",(mean+std), "diff1",1.0 - (mean+(sigma*std)), "diff2",diff)

    except:
        #print(sys.exc_info())
        #traceback.print_exc(file=sys.stdout)
        print("A", a)
        print("B", b)
        print("a-b", abs(a - b))
        print("maxone", max_one)
        print("Z",Z)
        print("diff",diff)

    return diff

def t_student(a,b,weights=None,p_value=0.05):
    aa = a.copy()
    bb = b.copy()
    if weights is not None:
        aa = (aa*weights).flatten()
        bb = (bb*weights).flatten()
    else:
        aa = aa.flatten()
        bb = bb.flatten()
    t2, p2 = scipy.stats.ttest_ind(aa, bb, equal_var=False)

    # print("t = " + str(t2))
    # print("p = " + str(p2))
    if p2<=p_value:
        #print("Significant: Means are statistically different. So they are not comparable")
        return False
    else:
        return True

def mann_whitney(a,b,weights=None, are_cosins=False, p_value=0.05, verbose=False):
    aa = a.copy()
    bb = b.copy()

    if not are_cosins:
        ddd = aa - bb
        max_one = a.copy()
        max_one[numpy.where(ddd < 0)] = bb[numpy.where(ddd < 0)]
        aa = aa / max_one
        bb = bb / max_one

    if weights is not None:
        aa = (aa*weights).flatten()
        bb = (bb*weights).flatten()
    else:
        aa = aa.flatten()
        bb = bb.flatten()
    try:
        t2, p2 = scipy.stats.mannwhitneyu(aa, bb)
    except ValueError:
        if verbose:
            print("Result is pvalue is more than threhsold")
            print()
            print(list(aa))
            print(list(bb))
            print("t = " + str(t2))
            print("p = " + str(p2))
            print()
        return True

    if p2<=p_value:
        # if verbose:
        #     print("Result is pvalue is less than threhsold")
        #     print()
        #     print(list(aa))
        #     print(list(bb))
        #     print("t = " + str(t2))
        #     print("p = " + str(p2))
        #     print()
        return False
    else:
        if verbose:
            print("Result is pvalue is more than threhsold")
            print()
            print(list(aa))
            print(list(bb))
            print("t = " + str(t2))
            print("p = " + str(p2))
            print()
        return True

def get_tm_score(reference, target):
    #### TM-score formula http://www.blopig.com/blog/2017/01/tm-score/
    assert(len(reference)==len(target))
    if not isinstance(reference[0],numpy.ndarray):
        reference = numpy.array([x.get_coord() for x in reference])
    if not isinstance(target[0], numpy.ndarray):
        target = numpy.array([x.get_coord() for x in target])

    diffm1 = [numpy.sqrt(numpy.sum((reference[t]-target[t])**2)) for t in range(len(reference))]

    d0 = 1.24 * numpy.cbrt(len(reference) - 15) - 1.8
    TM_Score1 = (1.0 / (len(reference))) * (numpy.sum([1.0 / (1.0 + (di / d0) ** 2) for di in diffm1]))
    return TM_Score1

def get_CA_distance_dictionary(pdb_model1, pdb_model2, max_rmsd=0.5, last_rmsd=1.0, recompute_rmsd=True, cycles=3, cycle=1, before_apply=None, get_superposed_atoms=False, force_reference_residues=False, data_tuple = None):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param pdb_model1:
    :param pdb_model2:
    :param max_rmsd:
    :param last_rmsd:
    :param recompute_rmsd:
    :param cycles:
    :param cycle:
    :param before_apply:
    :param get_superposed_atoms:
    :param force_reference_residues:
    :param data_tuple:
    :return:
    """
    if (isinstance(pdb_model1,str) and os.path.exists(pdb_model1)) or isinstance(pdb_model1,io.StringIO):
        if isinstance(pdb_model1,io.StringIO):
            pdb_model1.seek(0)
        structure1 = get_structure("a",pdb_model1)
        list_CA1 = [residue['CA'] for residue in Bio.PDB.Selection.unfold_entities(structure1,'R') if residue.has_id("CA")]
    else:
        list_CA1 = pdb_model1

    if (isinstance(pdb_model2,str) and os.path.exists(pdb_model2)) or isinstance(pdb_model2,io.StringIO):
        if isinstance(pdb_model2,io.StringIO):
            pdb_model2.seek(0)
        structure2 = get_structure("b",pdb_model2)
        list_CA2 = [residue['CA'] for residue in Bio.PDB.Selection.unfold_entities(structure2,'R') if residue.has_id("CA")]
    else:
        list_CA2 = pdb_model2

    if cycle == cycles:
        max_rmsd = last_rmsd
    #print("PDBMODEL1",type(pdb_model1),"PDBMODEL2",type(pdb_model2),"len(1)",len(list_CA1),"len(2)",len(list_CA2),"type(1)",type(list_CA1),"type(2)",type(list_CA2))
    rmsd = 100
    if before_apply is not None and isinstance(before_apply,tuple) or isinstance(before_apply,list):
        R, t = before_apply
        # print("---",full_ca[0].get_full_id(),full_ca[0].get_coord())
        #NOTE: it is essential that R,t corresponds to the transformation of the pdb_model2 onto pdb_mopdel1
        list_CA2 = transform_atoms(list_CA2, R, t)
    elif before_apply is not None and isinstance(before_apply,str) and before_apply == "automatic":
        allatoms_ana = [atm for res in list_CA2 for atm in res.get_parent()]
        list_CA2, rmsd, R, t = get_rmsd_and_RT(list_CA1, list_CA2, allatoms_ana)
        list_CA2 = [atom for atom in allatoms_ana if atom.get_name().lower() == "ca"]
        print("rmsd initial using all atoms is",rmsd)

    distance_hash = {ca1.get_full_id(): [(ca1-ca2,ca2.get_full_id(),ca1,ca2) for ca2 in list_CA2] for ca1 in list_CA1}
    if not force_reference_residues or data_tuple is None:
        distance_hash2 = {key:sorted(distance_hash[key], key=lambda x: x[0])[0] for key in distance_hash if (sorted(distance_hash[key], key=lambda x: x[0])[0][0] <= max_rmsd) or (sorted(distance_hash[key], key=lambda x: x[0])[0][0] <= max_rmsd+1.0 and
                                                                                                                                                                (key[0],key[1],key[2],(key[3][0],key[3][1]-1,key[3][2]),key[4]) in distance_hash and
                                                                                                                                                                              sorted(distance_hash[(key[0],key[1],key[2],(key[3][0],key[3][1]-1,key[3][2]),key[4])], key=lambda x: x[0])[0][0] <= max_rmsd and
                                                                                                                                                                              (key[0],key[1],key[2],(key[3][0],key[3][1]+1,key[3][2]),key[4]) in distance_hash and
                                                                                                                                                                              sorted(distance_hash[(key[0],key[1],key[2],(key[3][0],key[3][1]+1,key[3][2]),key[4])], key=lambda x: x[0])[0][0] <= max_rmsd)}
        todelete = []

        for i, key in enumerate(sorted(distance_hash2.keys())):
            prev = False
            post = False
            if i > 0 and (key[0], key[1], key[2], (key[3][0], key[3][1] - 1, key[3][2]), key[4]) in distance_hash2:
                prev = True
            if i < len(distance_hash2.keys()) - 1 and (
            key[0], key[1], key[2], (key[3][0], key[3][1] + 1, key[3][2]), key[4]) in distance_hash2:
                post = True

            if not (prev or post):
                todelete.append(key)

        for key in todelete:
            del distance_hash2[key]
    else:
        distance_hash2 = {key:sorted(distance_hash[key], key=lambda x: x[0]) for key in distance_hash if (sorted(distance_hash[key], key=lambda x: x[0])[0][0] <= max_rmsd) or (sorted(distance_hash[key], key=lambda x: x[0])[0][0] <= max_rmsd+1.0 and  (key[0],key[1],key[2],(key[3][0],key[3][1]-1,key[3][2]),key[4]) in distance_hash and
                                                                                                                                                                              sorted(distance_hash[(key[0],key[1],key[2],(key[3][0],key[3][1]-1,key[3][2]),key[4])], key=lambda x: x[0])[0][0] <= max_rmsd and
                                                                                                                                                                              (key[0],key[1],key[2],(key[3][0],key[3][1]+1,key[3][2]),key[4]) in distance_hash and
                                                                                                                                                                              sorted(distance_hash[(key[0],key[1],key[2],(key[3][0],key[3][1]+1,key[3][2]),key[4])], key=lambda x: x[0])[0][0] <= max_rmsd)}

        #1: Find the contact points between paired strands:
        (graph_ref, graph_targ, associations_list, map_reference, map_target) = data_tuple

        map_target = {k[2:4]:v for k,v in map_target.items()}
        map_reference = {k[2:4]:v for k,v in map_reference.items()}

        possibi = []
        for associations in associations_list:
            centers_ref = []
            corresp = []
            mins = []
            for fragr in graph_ref.vs:
                fragrl = [tuple(tr[2:4]) for tr in fragr["reslist"]]
                namer = fragr["name"]
                corre = associations[namer]
                lipol = []
                for key in distance_hash2.keys():
                    if key[2:4] in fragrl:
                        for q,value in enumerate(distance_hash2[key]):
                            if map_target[value[1][2:4]] == corre:
                                lipol.append((value,key[2:4]))
                                break
                lipo = sorted(lipol, key=lambda x:x[0][0])[0]
                diffe = [abs(lipol[c+1][0][0]-lipol[c][0][0]) for c in range(len(lipol)-1)]
                #print("DIFFE",diffe)
                where = fragrl.index(lipo[1])
                centers_ref.append(((0,where,len(fragrl)), fragr["reslist"][where]))
                corresp.append((centers_ref[-1],lipo[0]))
                mins.append(lipo[0][0])
            print("ASSO",associations)
            print("MINS",mins)
            possibi.append((sum(mins)/len(mins),associations,centers_ref,corresp))

        (score_possi, associations, centers_ref, corresp) = sorted(possibi, key=lambda x: x[0])[0]
        #centers_ref = [((0,round(int(len(frag["reslist"])/2)),len(frag["reslist"])),frag["reslist"][round(int(len(frag["reslist"])/2))]) for frag in graph_ref.vs]
        map_secstr_ref = {tuple(key[1:4]):frag.index for frag in graph_ref.vs for key in [tuple(k[1:4]) for k in frag["reslist"]]}
        map_secstr_targ = {tuple(key[1:4]):frag.index for frag in graph_targ.vs for key in [tuple(k[1:4]) for k in frag["reslist"]]}
        #corresp = [(center,distance_hash2[tr][0]) for center in centers_ref for tr in distance_hash2.keys() if tuple(tr[2:4])==tuple(center[1][2:4])]

        #print("CENTERS_REF",centers_ref)
        #print("MAP SECSTR REF",map_secstr_ref)
        #print("MAP SECSTR TARG",map_secstr_targ)

        toadd = []
        for corre in corresp:
            lista_res_targ = graph_targ.vs[map_secstr_targ[tuple(corre[1][1][2:4])]]["reslist"]
            central = [c for c,d in enumerate(lista_res_targ) if tuple(d[2:4])==tuple(corre[1][1][2:4])][0]
            sublist_ref = graph_ref.vs[map_secstr_ref[tuple(corre[0][1][2:4])]]["reslist"]
            start = central-corre[0][0][1]
            fine = central+(corre[0][0][2]-corre[0][0][1])
            add_at_start = 0
            add_at_fine = 0
            if start < 0:
                add_at_start = abs(start)
                start = 0
            if fine > len(lista_res_targ):
                add_at_fine = fine-len(lista_res_targ)
                fine = len(lista_res_targ)
            sublist_tar = lista_res_targ[start:fine]
            if len(sublist_ref) != len(sublist_tar):
                # print("CORRE", corre)
                # print("LISTA RES TARG", lista_res_targ)
                # print("CENTRAL", central)
                # print("SUBLIST REF", sublist_ref)
                # print("START", start, "FINE", fine)
                # print("START", start, "FINE", fine, "START ADD", add_at_start, "FINE ADD", add_at_fine)
                # print("SUBLIST TARG", sublist_tar)
                if add_at_start > 0 and add_at_fine > 0:
                    sublist_tar = sublist_tar+[sublist_tar[-1] for _ in range(add_at_start)]
                    sublist_tar = sublist_tar+[sublist_tar[-1] for _ in range(add_at_fine)]
                elif add_at_start == 0:
                    libre = start
                    if add_at_fine > libre:
                        sublist_tar = sublist_tar+lista_res_targ[0:libre]
                        sublist_tar = sublist_tar+[sublist_tar[-1]for _ in range(add_at_fine-libre)]
                    else:
                        sublist_tar = sublist_tar+lista_res_targ[libre-add_at_fine:libre]
                elif add_at_fine == 0:
                    libre = len(lista_res_targ)-fine
                    if add_at_start > libre:
                        sublist_tar = sublist_tar+lista_res_targ[fine:]
                        sublist_tar = sublist_tar+[sublist_tar[-1] for _ in range(add_at_start-libre)]
                    else:
                        sublist_tar = sublist_tar+lista_res_targ[fine:fine+add_at_start]

                #print("SUBLIST TAR CORRECTED",sublist_tar)
                if len(sublist_ref) != len(sublist_tar):
                    print("ERROR: they should be of the same size")
                    print(sublist_ref)
                    print(sublist_tar)
                    sys.exit(1)
            toadd.append(zip(sublist_ref,sublist_tar))

        for iterate in toadd:
            for sr,st in iterate:
                found = False
                for key in distance_hash2:
                    if key[1:4] == tuple(sr[1:4]):
                        for q,vlo in enumerate(distance_hash2[key]):
                            if vlo[1][1:4] == tuple(st[1:4]):
                                distance_hash2[key] = distance_hash2[key][q]
                                found = True
                                break
                        if found:
                            break

        # for k,v, in distance_hash2.items():
        #      print(k,v)
        # quit()

    distance_hash = distance_hash2
    # for key,value in distance_hash.items():
    #      print("--",key,value)

    allAtoms = []

    if recompute_rmsd:
        l1 = sorted(distance_hash.keys(), key=lambda x: x[1])
        lit1 = [distance_hash[c1][2] for c1 in l1]
        lit2 = [distance_hash[c1][3] for c1 in l1]
        if len(lit1) > 0 and len(lit1) == len(lit2):
            allatoms_ana = [atm for res in list_CA2 for atm in res.get_parent()]
            allAtoms, rmsd, R, t = get_rmsd_and_RT(lit1, lit2, allatoms_ana)
            allAtoms = [atom for atom in allatoms_ana if atom.get_name().lower() == "ca"]
            return get_CA_distance_dictionary(list_CA1, allAtoms, max_rmsd=max_rmsd if cycle+1 <= cycles else last_rmsd,
                                              recompute_rmsd=True if cycle+1 <= cycles else False, cycle=cycle+1,
                                              get_superposed_atoms=get_superposed_atoms)

    if get_superposed_atoms:
        l1 = sorted(distance_hash.keys(), key=lambda x: x[1])
        lit1 = [distance_hash[c1][2] for c1 in l1]
        lit2 = [distance_hash[c1][3] for c1 in l1]
        distance_hash = {key: distance_hash[key][:-2] for key in distance_hash}
        if len(lit1) > 0 and len(lit1) == len(lit2):
            allatoms_ana = [atm for res in list_CA2 for atm in res.get_parent()]
            allAtoms, rmsd, R, t = get_rmsd_and_RT(lit1, lit2, allatoms_ana)
            allAtoms = [atom for atom in allatoms_ana if atom.get_name().lower() == "ca"]
            return distance_hash, allAtoms, rmsd
        else:
            return distance_hash, list_CA2, rmsd
    else:
        distance_hash = {key: distance_hash[key][:-2] for key in distance_hash}
        return distance_hash

def get_rmsd_and_RT(ca_list1,ca_list2,full_ca, transform=True, n_iter=5):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param ca_list1: list of Biopython atom objects
    :type ca_list1: list
    :param ca_list2: list of Biopython atom objects
    :type ca_list2: list
    :param full_ca:
    :type full_ca:
    :return: allAtoms,rmsd,R,t in case transform is True, rmsd,R,t otherwise
    :return allAtoms:
    :rtype :
    :return rmsd:
    :rtype :
    :return R:
    :rtype :
    :return t:
    :rtype :
    """
    #structure1 = get_structure("cmp2", reference)
    #full_ca = [residue['CA'] for residue in Bio.PDB.Selection.unfold_entities(structure1, 'R')]

    atom_list_a = [atom.get_coord() if hasattr(atom,"get_coord") else atom for atom in ca_list1]
    atom_list_b = [atom.get_coord() if hasattr(atom,"get_coord") else atom for atom in ca_list2]
    atom_list_a = numpy.asarray(atom_list_a)
    atom_list_b = numpy.asarray(atom_list_b)
    transf, rmsd_list = csb.bio.utils.fit_wellordered(atom_list_a, atom_list_b, n_iter=n_iter, full_output=True,
                                                             n_stdv=2, tol_rmsd=0.005, tol_stdv=0.0005)
    if len(rmsd_list) > 0:
        rmsd = rmsd_list[-1][1]
    else:
        rmsd = 100

    R, t = transf
    #print("---",full_ca[0].get_full_id(),full_ca[0].get_coord())
    if transform:
        allAtoms = transform_atoms(full_ca, R, t)
        return allAtoms,rmsd,R,t
    else:
        return rmsd,R,t


def get_rmsd(coords1,coords2):
    """ Computes the RMSD between two equal sets of atoms

    :param coords1: numpy array that contains the biopython atom.get_coord() for each atom
    :type coords1: numpy array
    :param coords2:
    :type coords2: numpy array
    :return:
    """
    diff = coords1 - coords2
    return numpy.sqrt(numpy.sum(numpy.sum(diff * diff, axis=1)) / coords1.shape[0])


def get_rmsd_from_distance_hash(distance_hash):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param distance_hash:
    :type distance_hash:
    :return: rmsd
    :rtype: float
    """
    summatory = 0
    for key in distance_hash.keys():
        summatory=summatory+distance_hash[key][0]
    n = len(distance_hash)
    rmsd = numpy.sqrt(summatory/n)
    return rmsd

def percentage_of_secondary_structure(g, n_res):
    nres_ah = 0
    nres_bs = 0
    for vertex in g.vs:
        if vertex['sstype'] == 'ah':
            nres_ah += len(vertex['reslist'])
        elif vertex['sstype'] == 'bs':
            nres_bs += len(vertex['reslist'])

    perc_ah = int(round((nres_ah / n_res) * 100))
    perc_bs = int(round((nres_bs / n_res) * 100))
    print("Percentage of residues annotated as ah is:", perc_ah)
    print("Percentage of residues annotated as bs is:", perc_bs)
    return(perc_ah, perc_bs)

def generate_secondary_structure_record(g):
    """ Generate the HELIX/SHEET record in pdb format from a graph. First, will be described all the helices and then
    the beta sheets. Beta-strands of the same sheet will be sorted by their position in space, not by their sequence.

    :author: Ana del Rocío Medina Bernal
    :email: ambcri@ibmb.csic.es
    :param g: graph where each node is a secondary structure and the edges are the geometrical relations between them
    :type g: igraph.Graph
    :return: pdbsearchin
    :rtype pdbsearchin: str
    """

    pdbsearchin = ""
    ah = sorted((x for x in g.vs if x['sstype'] == 'ah'), key=lambda y: y['reslist'][0][3][1])

    try:
        bs_0 = sorted((x for x in g.vs if x['sstype'] == 'bs'),
                      key=lambda y: y['sheet'])  # Try to sort by number of sheet

    except: #If there is some sheet=None or sheet attribute does not exist
        bs_0 = [x for x in g.vs if x['sstype'] == 'bs']
        no_sheet= 1000
        none_sheet = 2000
        for fragment in bs_0:
            if not 'sheet' in fragment.attributes(): # If no sheet attribute
                fragment['sheet'] = no_sheet
                no_sheet += 1
            if fragment['sheet'] == None: # If sheet attribute is None
                fragment['sheet'] = none_sheet
                none_sheet += 1
    bs_0 = sorted((x for x in bs_0 if x['sstype'] == 'bs'),
                  key=lambda y: y['sheet'])  # Sort in case there is a bs with a fragment sheet different from 1000

    coil = sorted((x for x in g.vs if x['sstype'] == 'coil'), key=lambda y: y['reslist'][0][3][1])

    if len(bs_0) > 0:
        current_sheet = bs_0[0]['sheet']
    bs = []
    bs_current = []

    for frag in bs_0:  # Sort by sheet id
        if frag['sheet'] == current_sheet:
            bs_current.append(frag)
        else:
            bs.append(bs_current)
            bs_current = [frag]
            current_sheet = frag['sheet']

    bs.append(bs_current) if len(bs_current) > 0 else print('There are no beta-strands')

    final_bs = []

    if len(bs) > 0:
        for sheet in bs:
            n_con = []
            for frag in sheet:  # Calculate number of neighbours in a radious of 10A
                vecino = [x for x in (g.neighborhood(frag.index)) if (
                    g.vs(x)['sheet'] == g.vs(frag.index)['sheet'] and x != frag.index and
                    g.es[g.get_eid(frag.index, x)]['mean'][1] < 10)] #NOTE: Another version of the code had 10 here
                n_con.append([frag, len(vecino)])
            n_con = sorted(n_con, key=lambda y: (y[1], y[0]['reslist'][0][3][1]))

            bs_sort = [[n_con[0][0], 0]]
            n_con.pop(0)

            while len(n_con) > 0:
                for i, element in enumerate(n_con):
                    try:
                        n_con[i][2] = (g.es[g.get_eid(bs_sort[-1][0], element[0])]['mean'])
                    except:
                        n_con[i].append(g.es[g.get_eid(bs_sort[-1][0], element[0])]['mean'])
                    try:
                        if ALEPH.BS_UD_EA[int(g.es[g.get_eid(bs_sort[-1][0], element[0])]['mean'][0])] > ALEPH.BS_UU_EA[
                            int(g.es[g.get_eid(bs_sort[-1][0], element[0])]['mean'][0])]:
                            n_con[i][3] = -1
                        else:
                            n_con[i][3] = 1
                    except:
                        if ALEPH.BS_UD_EA[int(g.es[g.get_eid(bs_sort[-1][0], element[0])]['mean'][0])] > ALEPH.BS_UU_EA[
                            int(g.es[g.get_eid(bs_sort[-1][0], element[0])]['mean'][0])]:
                            n_con[i].append(-1)
                        else:
                            n_con[i].append(1)

                n_con = sorted(n_con, key=lambda y: y[2][1])
                bs_sort.append([n_con[0][0], n_con[0][3]])
                n_con.pop(0)
            final_bs.append(bs_sort)

    serNum = 1
    for frag in ah:
        pdbsearchin += get_helix_line(serNum, "H", ATOAAAMAP[frag["sequence"][0]],
                                      frag["reslist"][0][2], frag["reslist"][0][3][1], frag["reslist"][0][3][2],
                                      ATOAAAMAP[frag["sequence"][-1]],
                                      frag["reslist"][-1][2], frag["reslist"][-1][3][1],
                                      frag["reslist"][-1][3][2], 1, "", len(frag["sequence"]))
        serNum += 1

    for sheet in final_bs:
        strand = 1
        for frag in sheet:
            pdbsearchin += get_sheet_line(strand, "B", len(sheet), ATOAAAMAP[frag[0]["sequence"][0]],
                                          frag[0]["reslist"][0][2], frag[0]["reslist"][0][3][1],
                                          frag[0]["reslist"][0][3][2], ATOAAAMAP[frag[0]["sequence"][-1]],
                                          frag[0]["reslist"][-1][2], frag[0]["reslist"][-1][3][1],
                                          frag[0]["reslist"][-1][3][2], frag[1])
            strand += 1

    return pdbsearchin

def get_structure(name, pdbf, get_header=False):
    """
    Parse and generate a Structure object from a PDB file.

    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param name: name id of the structure
    :type name: str
    :param pdbf: path of the pdb file
    :type pdbf: Union[str,io.TextIOWrapper]
    :return: structure object
    :rtype: Bio.PDB.Structure
    """
    if isinstance(pdbf, io.StringIO):
        pdbf.seek(0)
    elif not os.path.exists(pdbf):
        pdbf = io.StringIO(SystemUtility.py2_3_unicode(pdbf))

    parser = Bio.PDB.PDBParser()
    structure = parser.get_structure(name, pdbf)

    if not get_header:
        return structure

    remarks = ""
    if isinstance(pdbf, io.StringIO):
        pdbf.seek(0)
        remarks = "".join([line for line in pdbf.read().splitlines() if line.startswith("REMARK")])
    elif not os.path.exists(pdbf):
        pdbf = io.StringIO(SystemUtility.py2_3_unicode(pdbf))
        remarks = "".join([line for line in pdbf.splitlines() if line.startswith("REMARK")])
    else:
        with open(pdbf,"r") as c:
            remarks = "".join([line for line in c.read().splitlines() if line.startswith("REMARK")])

    return structure,remarks

def rename_hetatm_and_icode(pdbf):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param pdbf:
    :return:
    """
    allreadlines = None
    with open(pdbf,"r") as f:
        allreadlines = f.readlines()

    for i,line in enumerate(allreadlines):
        if line.startswith("HETATM"):
            line = "ATOM  "+line[6:]
        if line.startswith("HETATM") or line.startswith("ATOM"):
            line = line[:26]+" "+line[27:]
        allreadlines[i] = line

    pdb = "".join(allreadlines)

    with open(pdbf,"w") as f:
        f.write(pdb)

def get_residue(structure, model, chain, residue):
    """
    Retrieve a specific residue from a structure
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure: structure object
    :type structure: Bio.PDB.Structure
    :param model: number of the model id
    :type model: int
    :param chain: chain name
    :type chain: chr
    :param residue: residue id
    :type residue: tuple
    :return: residue object if found or None
    :rtype: Bio.PDB.Residue or None
    """
    if isinstance(residue, tuple):
        if model is not None:
            lir = [resil for resil in Bio.PDB.Selection.unfold_entities(structure[model][chain], "R") if resil.get_id() == residue]
        else:
            lir = [resil for m in structure for resil in Bio.PDB.Selection.unfold_entities(m[chain], "R") if resil.get_id() == residue]
    else:
        if model is not None:
            lir = [resil for resil in Bio.PDB.Selection.unfold_entities(structure[model][chain], "R") if
                   resil.get_id()[1] == residue]
        else:
            lir = [resil for m in structure for resil in Bio.PDB.Selection.unfold_entities(m[chain], "R") if
                   resil.get_id()[1] == residue]

    if len(lir) > 0:
        return lir[0]
    else:
        # print(model,chain,residue)
        # print()
        # for resil in Bio.PDB.Selection.unfold_entities(structure[model][chain], "R"):
        #     print(resil.get_full_id())
        return None

def get_dictio_resi_to_secstr(graph,structure):
    dic_resi = {tuple(resi[1:-1]):(frag["sstype"],frag.index) for frag in graph.vs for resi in frag["reslist"] if frag["sstype"] in ["ah","bs"]}

    for resi in get_list_of_residues(structure):
        if resi.has_id("CA") and resi.has_id("C") and resi.has_id("O") and resi.has_id("N") and tuple(resi.get_full_id()[1:]) not in dic_resi:
            dic_resi[tuple(resi.get_full_id()[1:])] = ("cc",None)

    return dic_resi

def get_list_of_residues(structure, model=None, sorting=False, interface = False):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure:
    :param model:
    :return:
    """
    lir = [resil for resil in Bio.PDB.Selection.unfold_entities(structure if model is None else structure[model], "R")]

    if sorting:
        lir = sorted(lir, key=lambda x: x.get_full_id()[3][2])
    if interface:
        lir = sorted(lir, key=lambda x: (x.get_full_id()[2], x.get_full_id()[3][2]))
    return lir

def get_list_of_atoms(structure, model=None, sorting=False, type_atom=None):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure:
    :type structure: Biopython structure object
    :param model:
    :type model:
    :param sorting: whether to sort or not the atoms list
    :type sorting: bool
    :param type_atom: whether to select particular atoms such as CA
    :type type_atom: str
    :return lir: list with the biopython atom objects
    :return type: list
    """
    lir = [resil for resil in Bio.PDB.Selection.unfold_entities(structure if model is None else structure[model], "A") if type_atom is None or resil.get_name()==type_atom.upper()]
    if sorting:
        lir = sorted(lir, key=lambda x: x.get_full_id())

    return lir

def get_atom(structure, model, chain, residue, atom):
    """
    Retrieve a specific residue from a structure
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure: structure object
    :type structure: Bio.PDB.Structure
    :param model: number of the model id
    :type model: int
    :param chain: chain name
    :type chain: chr
    :param residue: residue id
    :type residue: tuple
    :return: residue object if found or None
    :rtype: Bio.PDB.Residue or None
    """
    lir = [atoml for atoml in Bio.PDB.Selection.unfold_entities(structure[model][chain][residue], "A") if atoml.get_id() == atom]
    if len(lir) > 0:
        return lir[0]
    else:
        return None

def get_backbone(residue, without_CB=True):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param residue:
    :param without_CB:
    :return:
    """
    if without_CB or not residue.has_id("CB"):
        return [residue["CA"],residue["C"],residue["O"],residue["N"]]
    else:
        return [residue["CA"],residue["C"],residue["O"],residue["N"],residue["CB"]]

def fetch_pdb(pdbid, outputfile):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param pdbid:
    :param outputfile:
    :return:
    """
    tries = 0
    while 1:
        try:
            baseurl = "ftp://ftp.wwpdb.org/pub/pdb/data/structures/all/pdb/pdb" + pdbid + ".ent.gz"
            try:
                urllib.request.urlretrieve(baseurl, outputfile)
            except:
                urllib2.urlopen(baseurl, outputfile)
            finally:
                break
        except:
            SystemUtility.warning("Cannot download " + str(pdbid) + ".pdb . Try again...")
            tries += 1
            if tries > 10:
                break

    pdbf = outputfile
    if tries > 10:
        SystemUtility.warning("Cannot download " + str(pdbid) + ".pdb . Stop trying. Skipping...")
        # raise Exception("Cannot download "+str(pdb)+".pdb . Stop trying. Skipping...")
        # NOTE: TEMPORANEO PDB is not working today
        # pdbf = "/localdata/PDBDB_20160412/"+pdb[1:3]+"/pdb"+pdb+".ent.gz"

    if pdbf.endswith(".gz"):
        fileObj = gzip.GzipFile(pdbf, 'rb');
        fileContent = fileObj.read()
        fileObj.close()
        os.remove(pdbf)
        pdbf = pdbf[:-3]  # elimino estensione .gz
        fou = open(pdbf, "w")
        fou.write(fileContent)
        fou.close()
    if pdbf.endswith(".ent"):
        pdbf2 = pdbf[:-4]  # elimino estensione .ent
        pdbf2 = pdbf2 + ".pdb"
        os.rename(pdbf, pdbf2)
        pdbf = pdbf2
    if os.path.basename(pdbf).startswith("pdb"):
        root, fileu = os.path.split(pdbf)
        pdbf2 = os.path.basename(pdbf)[3:]  # elimino la parola pdb
        pdbf2 = os.path.join(root, pdbf2)
        os.rename(pdbf, pdbf2)
        pdbf = pdbf2
    shutil.move(pdbf, outputfile)

def write_pdb(structure,path):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure:
    :param path:
    :return:
    """
    io = Bio.PDB.PDBIO()
    io.set_structure(structure)
    io.save(path)

def distance_sq(X, Y):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param X:
    :param Y:
    :return:
    """
    return ((numpy.asarray(X) - Y) ** 2).sum(-1)

def transform_atoms(atoms, R, t, s=None, invert=False):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param atoms:
    :param R:
    :param t:
    :param s:
    :param invert:
    :return:
    """
    Y = []

    for atom in atoms:
        Y.append(atom.get_coord())

    Y = transform(Y, R, t, s=s, invert=invert)

    for i in range(len(Y)):
        y = Y[i]
        atoms[i].set_coord(y)

    return atoms

def transform(Y, R, t, s=None, invert=False):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param Y:
    :param R:
    :param t:
    :param s:
    :param invert:
    :return:
    """
    if invert:
        x = numpy.dot(Y - t, R)
        if s is not None:
            s = 1. / s
    else:
        x = numpy.dot(Y, R.T) + t
    if s is not None:
        x *= s
    return x

def check_continuity(res1, res2, swap=True, verbose=False):
    """
    Check if two residues are continous in a structure. It checks the 3d cordinates not the residue ids.
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param res1: First Residue
    :type res1: Bio.PDB.Residue
    :param res2: Following Residue
    :type res2: Bio.PDB.Residue
    :return: True if they are cntinous, False if they are not
    :rtype: bool
    """

    try:
        resaN = res2["N"]
        prevResC = res1["C"]
    except:
        return False

    checkCont = numpy.sqrt(((resaN.get_coord()-prevResC.get_coord())**2).sum(axis=0))

    if verbose:
        print("CHECKING",resaN.get_full_id(),prevResC.get_full_id(),checkCont)

    if checkCont <= 1.5:
        return True
    else:
        if swap:
            return check_continuity(res2, res1, swap=False)
        else:
            return False

def get_helix_line(idnumber, idname, restartname, chainstart, restartnumber, resicodestart, resendname, chainend, resendnumber, resicodeend, typehelix, comment, lenhelix):
    """
    For an helix it returns the string formatted representing the HELIX record in the pdb
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param idnumber: index number
    :type idnumber: int
    :param idname: idname
    :type idname: str
    :param restartname: first residue name
    :type restartname: str
    :param chainstart: chain of the first residue
    :type chainstart: chr
    :param restartnumber: number of the first residue
    :type restartnumber: int
    :param resicodestart: icode of the first residue
    :type resicodestart: str
    :param resendname:  name of the last residue
    :type resendname: str
    :param chainend:  chain of the last residue
    :type chainend: chr
    :param resendnumber:  number of the last residue
    :type resendnumber: int
    :param resicodeend: icode of the last residue
    :type resicodeend: str
    :param typehelix:  type of helix (see PDB standard for HELIX)
    :type typehelix: int
    :param comment:  string comment
    :type comment: str
    :param lenhelix: number of reasidues in the helix
    :type lenhelix: int
    :return record: HELIX record
    :rtype record: str
    """

    # HELIX_FORMAT_STRING = "{:<6s} {:>3d} {:>3s} "
    # print('llega', idnumber,idname,restartname,chainstart,restartnumber,resicodestart,resendname,chainend,resendnumber,resicodeend,typehelix,comment,lenhelix)
    HELIX_FORMAT_STRING = "{:<6s} {:>3d} {:>3s} {:>3s} {:1s} {:>4d}{:1s} {:>3s} {:1s} {:>4d}{:1s}{:>2d}{:>30s} {:>5d}"
    return HELIX_FORMAT_STRING.format("HELIX",idnumber,idname,restartname,chainstart,restartnumber,resicodestart,resendname,chainend,resendnumber,resicodeend,typehelix,comment,lenhelix)+"\n"


def get_sheet_line(idnumber, idnamesheet, nofstrandsinsheet, restartname, chainstart, restartnumber, resicodestart, resendname, chainend, resendnumber, resicodeend, sense):
    """
    For a beta strand it returns the string formatted representing the SHEET record in the pdb
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param idnumber: index number
    :type idnumber: int
    :param idnamesheet: idname
    :type idnamesheet: str
    :param nofstrandsinsheet: number of beta strand in the beta sheet
    :type nofstrandsinsheet: int
    :param restartname: first residue name
    :type restartname: str
    :param chainstart: chain of the first residue
    :type chainstart:  chr
    :param restartnumber: number of the first residue
    :type restartnumber: int
    :param resicodestart: icode of the first residue
    :type resicodestart: str
    :param resendname: name of the last residue
    :type resendname: str
    :param chainend: chain for the last residue
    :type chainend: chr
    :param resendnumber: number of the last residue
    :type resendnumber: int
    :param resicodeend:  icode of the last residue
    :type resicodeend: str
    :param sense: direction of the strand (see PDB standard for SHEET)
    :type sense: int
    :return record:  SHEET record
    :rtype record: str
    """

    SHEET_FORMAT_STRING = "{:<6s} {:>3d} {:>3s}{:>2d} {:>3s} {:1s}{:>4d}{:1s} {:>3s} {:1s}{:>4d}{:1s}{:>2d} {:>4s}{:>3s} {:1s}{:>4s}{:1s} {:>4s}{:>3s} {:1s}{:>4s}{:1s}"
    return SHEET_FORMAT_STRING.format("SHEET",idnumber,idnamesheet,nofstrandsinsheet,restartname,chainstart,restartnumber,resicodestart,resendname,chainend,resendnumber,resicodeend,sense,"","","","","","","","","","")+"\n"

def get_atom_line(atom, element, hetfield, segid, atom_number, resname, resseq, icode, chain_id, normalize=False, bfactorNor=25.00, bfactorNorSidechains=50.00, charge=" ", applyRt=None, predicted_model=False):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param atom:
    :param element:
    :param hetfield:
    :param segid:
    :param atom_number:
    :param resname:
    :param resseq:
    :param icode:
    :param chain_id:
    :param normalize:
    :param bfactorNor:
    :param charge:
    :return:
    """

    ATOM_FORMAT_STRING = "%s%6i %-4s%c%3s %c%4i%c   %8.3f%8.3f%8.3f%6.2f%6.2f      %4s%2s%2s\n"

    if hetfield != " ":
        record_type = "HETATM"
    else:
        record_type = "ATOM "

    element = element.strip().upper()
    # print "ELEMENT value was: ",type(element),element
    element = element[0]
    # print "ELEMENT value is: ",type(element),element

    name = atom.get_fullname()
    altloc = atom.get_altloc()

    if applyRt is not None:
        coo = atom.get_coord()
        coo = numpy.array([coo])
        coo = transform(coo, applyRt[0], applyRt[1])
        x, y, z = coo[0]
    else:
        x, y, z = atom.get_coord()

    occupancy = atom.get_occupancy()

    if not normalize:
        bfactor = atom.get_bfactor()

    if normalize and not predicted_model:
        bfactor = bfactorNor

    if normalize and predicted_model:
        if atom.get_name().lower() in ["ca", "c", "o", "n", "cb"]:
            bfactor = bfactorNor
        else:
            bfactor = bfactorNorSidechains

    args = (record_type, atom_number, name, altloc, resname, chain_id, resseq, icode, x, y, z, occupancy, bfactor, segid, element, charge)
    ala = ATOM_FORMAT_STRING % args
    if record_type == "HETATM":
        clu = ala.split()
        spaziBianchi = 5 - len(clu[1])
        stri = "HETATM"
        for ulu in range(spaziBianchi):
            stri += " "
        stri += clu[1]
        ala = stri + ala[12:]
    return ala

def get_pdb_from_structure(structure,model):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param structure:
    :param model:
    :return:
    """
    atoms = get_list_of_atoms(structure,model)
    return get_pdb_from_list_of_atoms(atoms, renumber=False, uniqueChain=False, chainId="A", chainFragment=False, diffchain=None, polyala=True, maintainCys=False, normalize=False, sort_reference=True)[0]

#NOTE: When chains with more than one character will be provided this method will rename everything to chain A
def get_pdb_from_list_of_atoms(reference, renumber=False, uniqueChain=False, chainId="A", chainFragment=False,
                               diffchain=None, polyala=True, maintainCys=False, normalize=False, sort_reference=True,
                               remove_non_res_hetatm=False, dictio_chains={}, write_pdb=False, path_output_pdb='',
			       applyRt=None, interface=False, predicted_model=False):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param reference:
    :param renumber:
    :param uniqueChain:
    :param chainId:
    :param chainFragment:
    :param diffchain:
    :param polyala:
    :param maintainCys:
    :param normalize:
    :param sort_reference:
    NOTE CM: these parameters have been added by me
    :param remove_non_res_hetatm:
    :param dictio_chains:
    :param write_pdb:
    :return:
    """
    pdbString = ""
    numero = 1
    resn = {}
    nur = 1
    lastRes = None
    prevChain = 0
    lich = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "Z",
            "J", "K", "X", "Y", "W", "a", "b", "c", "d", "e", "f", "g", "h", "i", "l", "m", "n", "o", "p", "q", "r",
            "s", "t", "u", "v", "z", "j", "k", "x", "y", "w", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
    # if chainFragment or uniqueChain:
    #   renumber = True

    if remove_non_res_hetatm:
        reference = [atm for atm in reference if atm.get_parent().get_resname().upper() in AAList]

    if not polyala:
        reference = [atm for res in reference for atm in res.get_parent()]
    elif maintainCys:
        sort_reference = True
        ref2 = []
        # reference = [atm if res.get_resname().lower() == "cys" else res for res in reference for atm in res.get_parent()]
        for i in range(len(reference)):
            atm = reference[i]
            if atm.get_name().lower() != "ca":
                continue
            if atm.get_parent().get_resname().lower() == "cys":
                for atm2 in atm.get_parent():
                    ref2.append(atm2)
            else:
                for atm2 in atm.get_parent():
                    if atm2.get_name().lower() in ["ca", "c", "o", "n", "cb"]:
                        ref2.append(atm2)
        reference = ref2

    elif polyala:
        reference = [atm for atm in reference if atm.get_name().lower() in ["ca", "c", "o", "n", "cb"]]

    if predicted_model:
        reference = [atm for atm in reference if not atm.get_name().lower().startswith("h") if not atm.get_name()[0].isdigit()]

    if sort_reference:
        reference = Bio.PDB.Selection.uniqueify(reference)
        # listore = sorted(reference, key=lambda x: x.get_full_id()) Before insertion code
        listore = sorted(reference, key=lambda x: (x.get_full_id()[1], x.get_full_id()[2], x.get_full_id()[3][2], x.get_full_id()[3][1], x.get_full_id()[4][0]))
    else:
        listore = reference

    if interface:
        for item in listore:
            hetfield, resseq, icode = item.get_parent().get_id()
            segid = item.get_parent().get_segid()
            chain_id = item.get_parent().get_parent().get_id()
            hetfield = " "

            if len(chain_id) > 1:
                uniqueChain = True
                chainId = "A"

            if item.get_parent().get_full_id() not in resn.keys():
                if lastRes != None:
                    if not check_continuity(lastRes.get_parent(), item.get_parent()):
                        if renumber:
                            # print("LASTRES", lastRes.get_parent().get_full_id(), "ITEM:", item.get_parent().get_full_id())
                            nur += 10
                        if chainFragment:
                            prevChain += 1
                new_chain_id = chain_id
                if dictio_chains != {}:
                    new_chain_id = dictio_chains[item.get_parent().get_full_id()]
                    # print('ew', item.get_full_id(), new_chain_id)

                elif uniqueChain:
                    new_chain_id = chainId
                elif chainFragment:
                    new_chain_id = lich[prevChain]
                if renumber:
                    resn[item.get_parent().get_full_id()] = (nur, new_chain_id)
                    # resn[(resseq, chain_id)] = (nur, new_chain_id)
                else:
                    resn[item.get_parent().get_full_id()] = (resseq, new_chain_id)
                    # resn[(resseq, chain_id)] = (resseq, new_chain_id)
                lastRes = item
                nur += 1
            tuplr = resn[item.get_parent().get_full_id()]
            resseq = tuplr[0]
            chain_id = tuplr[1]
            icode = " "
            orig_atom_num = numero
            numero += 1

            resname = item.get_parent().get_resname()
            element = item.get_name()
            pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id,
                                       normalize=normalize, applyRt=applyRt, predicted_model=predicted_model)

    else:
        for item in listore:

            #if item.get_full_id() in set(seen):
            #    continue

            #seen.append(item.get_full_id())
            orig_atom_num = item.get_serial_number()
            hetfield, resseq, icode = item.get_parent().get_id()
            segid = item.get_parent().get_segid()
            chain_id = item.get_parent().get_parent().get_id()
            hetfield = " "

            # NOTE: PDB cannot support chains with more than one character in that case a pdb could not be written.
            # For now I am just renaming to chain A, but possibly we should move to .cif
            if len(chain_id) > 1:
                uniqueChain = True
                chainId = "A"

            if (resseq, chain_id) not in resn.keys():
                if lastRes != None:
                    # print "Checking Continuity",lastRes.get_parent().get_full_id(),item.get_parent().get_full_id(),checkContinuity(lastRes.get_parent(),item.get_parent())
                    # print "Checking Continuity",item.get_parent().get_full_id(),lastRes.get_parent().get_full_id(),checkContinuity(item.get_parent(),lastRes.get_parent())
                    # print lich[prevChain]
                    # print
                    if not check_continuity(lastRes.get_parent(), item.get_parent()):
                        if renumber:
                            #print("LASTRES", lastRes.get_parent().get_full_id(), "ITEM:", item.get_parent().get_full_id())
                            nur += 10
                        if chainFragment:
                            prevChain += 1
                new_chain_id = chain_id
                if dictio_chains!={}:
                    new_chain_id=dictio_chains[item.get_parent().get_full_id()]
                    # print('ew', item.get_full_id(), new_chain_id)

                elif uniqueChain:
                    new_chain_id = chainId
                elif chainFragment:
                    new_chain_id = lich[prevChain]
                if renumber:
                    resn[(resseq, chain_id)] = (nur, new_chain_id)
                else:
                    resn[(resseq, chain_id)] = (resseq, new_chain_id)
                lastRes = item
                nur += 1
            tuplr = resn[(resseq, chain_id)]
            resseq = tuplr[0]
            chain_id = tuplr[1]
            icode = " "
            orig_atom_num = numero
            numero += 1

            resname = item.get_parent().get_resname()
            element = item.get_name()
            pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id, normalize=normalize, applyRt=applyRt, predicted_model=predicted_model)

    if diffchain != None and len(diffchain) > 0:
        prevChain += 1
        lastRes = None
        for item in diffchain:
            orig_atom_num = item.get_serial_number()
            hetfield, resseq, icode = item.get_parent().get_id()
            segid = item.get_parent().get_segid()
            chain_id = item.get_parent().get_parent().get_id()
            hetfield = " "
            if (resseq, chain_id) not in resn.keys():
                if lastRes != None:
                    if not check_continuity(lastRes.get_parent(), item.get_parent()):
                        if renumber:
                            #print("LASTRES",lastRes.get_parent().get_full_id(),"ITEM:",item.get_parent().get_full_id())
                            nur += 10
                        if chainFragment:
                            prevChain += 1
                new_chain_id = chain_id
                if uniqueChain:
                    new_chain_id = chainId
                elif chainFragment:
                    new_chain_id = lich[prevChain]

                if renumber:
                    resn[(resseq, chain_id)] = (nur, new_chain_id)
                else:
                    resn[(resseq, chain_id)] = (resseq, new_chain_id)
                lastRes = item
                nur += 1
            tuplr = resn[(resseq, chain_id)]
            resseq = tuplr[0]
            chain_id = tuplr[1]
            icode = " "
            orig_atom_num = numero
            numero += 1

            resname = item.get_parent().get_resname()
            element = item.get_name()
            pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id, normalize=normalize, applyRt=applyRt, predicted_model=predicted_model)

    if write_pdb:
        fichmodel=open(path_output_pdb,'w')
        fichmodel.write(pdbString)
        fichmodel.close()
    return pdbString, resn


def get_pdb_from_list_of_frags(nModel, allfrags, structure, pathBase, dizioConv={}, externalRes=[], normalize=False):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param nModel:
    :param allfrags:
    :param structure:
    :param pathBase:
    :param dizioConv:
    :param externalRes:
    :param normalize:
    :return:
    """
    atoms = []
    for frag in allfrags:
        for resi in frag["reslist"]:
            atoms += get_backbone(get_residue(structure, resi[1], resi[2], resi[3]),without_CB=False)

    pdb_string, dict_conv = get_pdb_from_list_of_atoms(atoms,normalize=normalize)

    pdbid = str(structure.get_id())
    nomeFilefine = pathBase + "/" + pdbid + "_"

    for y in range(len(allfrags)):
        fragment = allfrags[y]
        for i in range(len(fragment["resIdList"])):
            re = fragment["resIdList"][i]
            if i == 0:
                if len(dizioConv.keys()) > 0:
                    nomeFilefine += str((dizioConv[(fragment["chain"], re, "CA")])[1]) + "*" + str(fragment["fragLength"])
                else:
                    nomeFilefine += str(re[1]) + str(re[2])
            if y != (len(allfrags) - 1) and i == 0:
                nomeFilefine += "_"
            elif y == (len(allfrags) - 1) and i == 0:
                if nModel != "":
                    nomeFilefine += "_" + nModel

    pdbString = ""
    pdbString += "REMARK TITLE " + nomeFilefine + "\n"
    pdb_string = pdbString+pdb_string

    return (nomeFilefine, pdb_string, dict_conv)


def read_pdb_ss_information_return_list_dics (input_pdb_file):
    """ Reads the pdb and collects the information about each HELIX and SHEET record in a list

    :author: Ana Medina
    :email: ambcri@ibmb.csic.es

    :param input_file_pdb: pdb file with the ss annotation
    :type input_pdb_file: str
    :return list_file_all: list containing a dictionary for each aa presented in a ss
    :rtype list_file_all: list
    :return False: In case no remark was found
    :rtype: bool
    """

    pdbfile = open(input_pdb_file)
    file_lines=pdbfile.readlines()
    pdbfile.close()
    list_file_all=[] #List with dictionaries. Each dictionary correspond to each ss element of the pdb file

    def catch_substring(a, b, line):
        try:
            return line[a:b].strip()
        except:
            return ""

    for line in file_lines:
        dic_var = {}
        dic_var['ss']= line[0:6].strip()
        if dic_var['ss'] == 'HELIX':
            dic_var['ser_num'] = int(catch_substring(7, 10, line)) #Serial number of the helix. This starts at 1 abd increases incrementally
            dic_var['helix_id'] = catch_substring(11, 14, line) #Helix identifier. In additio to a serial number, each helix is given and alphanumeric character helix identifier
            dic_var['init_res_name'] = catch_substring(15, 18, line) #Name of the initial residue
            dic_var['init_chain_id'] = catch_substring(19, 20, line) #Chain identifier for the chain containing this helix
            dic_var['init_seq_num'] = int(catch_substring(21, 25, line)) #Sequence number of the initial residue
            dic_var['init_i_code'] = catch_substring(25, 26, line) #Insertion code of the initial residue
            dic_var['end_res_name'] = catch_substring(27, 30, line) #Name of the terminal residue of the helix
            dic_var['end_chain_id'] = catch_substring(31, 32, line) #Chain identifier fot the chain containing this helix
            dic_var['end_seq_num'] = int(catch_substring(33, 37, line)) #Sequence number of the terminal residue
            dic_var['end_i_code'] = catch_substring(37, 38, line) #Insertion code of the terminal residue
            dic_var['helix_class'] = int(catch_substring(38, 40, line)) #Helix class (see ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_A4.pdf)
            dic_var['comment'] = catch_substring(40, 70, line) #Comment about this helix
            dic_var['length'] = int(catch_substring(71, 76, line)) #Length of this helix
            list_file_all.append(dic_var)
        elif dic_var['ss'] == 'SHEET':
            dic_var['strand'] = int(catch_substring(7, 10, line)) #Strand number which starsts aty 1 fot each strand within a sheet and increases by one
            dic_var['sheet_id'] = catch_substring(11, 14, line) #Sheet identifier
            dic_var['num_strands'] = int(catch_substring(14, 16, line)) #Number of strands in sheet
            dic_var['init_res_name'] = catch_substring(17, 20, line) #Residue name of initial residue
            dic_var['init_chain_id'] = catch_substring(21, 22, line) #Chain identifier of initial residue in strand
            dic_var['init_seq_num'] = int(catch_substring(22, 26, line)) #Sequence number of initial residue in strand
            dic_var['init_i_code'] = catch_substring(26, 27, line) #Insertion code of initial residue in strand
            dic_var['end_res_name'] = catch_substring(28, 31, line) #Residue name of terminal residue
            dic_var['end_chain_id'] = catch_substring(32, 33, line) #Chain identifier of terminal residue
            dic_var['end_seq_num'] = int(catch_substring(33, 37, line)) #Sequence number of terminal residue
            dic_var['end_i_code'] = catch_substring(37, 38, line) #Insertion code of terminal residue
            dic_var['sense'] = int(catch_substring(38, 40, line)) #Sense of strand with respect to previous strand in the sheet. 0 if first strand, 1 if parallel, and -1 if anti-parallel
            dic_var['cur_atom'] = catch_substring(41, 45, line) #Registration. Atom name in current strand
            dic_var['cur_res_name'] = catch_substring(45, 48, line) #Registration. Residue name in current strand
            dic_var['cur_chain_id'] = catch_substring(49, 50, line) #Registration. Chain identifier in current strand
            dic_var['cur_res_seq'] = catch_substring(50, 54, line) #Residue sequence number in current strand
            dic_var['cur_i_code'] = catch_substring(54, 55, line) #Registration. Insertion code in current strand
            dic_var['prev_atom'] = catch_substring(56, 60, line) #Registration. Atom name in previous strand
            dic_var['prev_res_name'] = catch_substring(60, 63, line) #Registration. Residue name in previous strand
            dic_var['prev_chain_id'] = catch_substring(64, 65, line) #Registration. Chain identifier in previous strand
            dic_var['prev_res_seq'] = catch_substring(65, 69, line) #Registration. Residue sequence number in previous strand
            dic_var['prev_i_code'] = catch_substring(69, 70, line) #Registration. Insertion code in previous strand
            list_file_all.append(dic_var)
        if len (list_file_all) == 0:
            return False
        else:
            return list_file_all


def reset_occupancies_in_pdb_to_one(pdb_filepath):
    """ Reads a pdb file of the disordered residues or atoms keeps the major occupancy one, setting them to one
    :author: Claudia Millan
    :email: cmncri@ibmb.csic.es
    :param pdb_filepath:
    :return: it just rewrites the pdb

    """
    oristru = get_structure(os.path.basename(pdb_filepath)[:-4],pdb_filepath)
    list_all_atoms = Bio.PDB.Selection.unfold_entities(oristru, 'A')
    list_atoms_to_write = []
    for atom in list_all_atoms:
        if atom.is_disordered():
            major_occ = 0
            for atomi in atom.disordered_get_list():
                #print('occupancy for atomi is ',atomi.occupancy)
                #print('id for atomi is ', atomi.get_full_id())
                if atomi.occupancy > major_occ:
                    atom_to_keep_id = atomi.get_full_id()
                    atom_to_keep = atomi
                    major_occ = atomi.occupancy
            atom_to_keep.occupancy = 1.0
            atom.set_altloc(' ')
            list_atoms_to_write.append(atom_to_keep)
        else:
            list_atoms_to_write.append(atom)

    #outpdbpath = pdb_filepath[:-4]+'_occreset.pdb' # only to check
    outpdbpath = pdb_filepath
    list_atoms_to_write = sorted(list_atoms_to_write,key=lambda x:x.get_parent().get_full_id()[3][1]) # sort by res number
    get_pdb_from_list_of_atoms(reference=list_atoms_to_write,renumber=False, uniqueChain=False, chainId="A",
                               chainFragment=False, polyala=False, maintainCys=False,
                               normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
                               write_pdb=True, path_output_pdb=outpdbpath)


def shredder_template_annotation(model_file, current_directory, bfacnorm=True, poliA=True, cys=False, remove_coil=True,
                                 nres_extend=0, min_alpha=7, min_beta=4, min_diff_ah=0.45, min_diff_bs=0.2,
                                 gyre_preserve_chains=False, algorithm_community='fastgreedy',
                                 pack_beta_community=False, homogenity_community=False, predicted_model=False):
    """ Annotates a model pdb file in terms of secondary and tertiary structure, using ALEPH.

    :author: Claudia Millan
    :email: cmncri@ibmb.csic.es


    :param model_file: path to the pdb file to process
    :type model_file: str
    :param current_directory: path of the current working directory
    :type current_directory: str
    :param bfacnorm: indicates whether to perform or not a bfactor normalization (default True)
    :type bfacnorm: boolean
    :param poliA: indicates whether to convert the model to only mainchain atoms (default True)
    :type poliA: boolean
    :param cys: indicates whether to leave or not cysteine residues untouched even if poliA is True (default False)
    :type cys: boolean
    :param remove_coil: indicates whether to leave the coil in the template or not (default True)
    :type remove_coil: boolean
    :param nres_extend: number of residues to add to secondary structure elements in the partial_coil case
    :type nres_extend: int
    :param min_alpha: minimum size in residues for any given helix in the template to be considered
    :type min_alpha: int
    :param min_beta: minimum size in residues for any given beta strand in the template to be considered
    :type min_beta: int
    :param min_diff_ah:
    :type min_diff_ah:
    :param min_diff_bs:
    :type min_diff_bs:
    :param gyre_preserve_chains:
    :type gyre_preserve_chains:
    :param algorithm_community:
    :type algorithm_community:
    :param pack_beta_community:
    :type pack_beta_community:
    :param homogenity_community:
    :type homogenity_community:

    :return dict_oristru: dictionary with keys being each residue number, and the following structure
             dict_template[nres] = {'residue_object': BioPython residue object,
                          'ori_full_id': tuple, 'ori_chain': str, 'ori_nres': int,'ori_nameres': str,
                          'ss_type_res': str (can be ah,bs,coil), 'ss_reslist': list, 'ss_id_res': str,
                          'first_ref_cycle_group': str, 'second_ref_cycle_group': str,'third_ref_cycle_group': str}
    :rtype dict_oristru: dict
    :return model_file: path to the modified template file
    :rtype model_file: str
    :return distance_matrix_CA:
    :rtype distance_matrix_CA: list of lists
    :return names_matrix:
    :rtypenames_matrix: list
    """

    print("\n * Info * ARCIMBOLDO_SHREDDER template treatment and annotation has started:")
    print("\n Processing ", model_file)

    # Prepare the input pdb file
    peptide_length = 3
    rename_hetatm_and_icode(model_file)
    reset_occupancies_in_pdb_to_one(model_file)
    parser = Bio.PDB.PDBParser()
    oristru_name = os.path.basename(model_file)[:-4]
    oristru = parser.get_structure(os.path.basename(model_file)[:-4], model_file)
    list_atoms = Bio.PDB.Selection.unfold_entities(oristru, 'A')

    dictio_chainid = {}
    # check the chain id for the case in which it might be empty
    set_chain_id = set([ atom.get_full_id()[2] for atom in list_atoms])
    if (len(set_chain_id)==1 and set_chain_id == set([' '])) or (len(set_chain_id) > 1 and not gyre_preserve_chains):
        for model in oristru:
            for chain in model:
                for residue in chain:
                    key_chain = residue.get_full_id()  # that is the key for the dictio_chains dictionary
                    try:
                        dictio_chainid[key_chain] = 'A'
                    except:
                        print('Some error happened')
                        sys.exit(0)

    if len(set_chain_id)==1 and not gyre_preserve_chains: # unique chain, we want to keep it the original
        unique_chain_id=(list(set_chain_id))[0]
        for model in oristru:
            for chain in model:
                for residue in chain:
                    key_chain = residue.get_full_id()  # that is the key for the dictio_chains dictionary
                    #print('SHERLOCK key_chain',key_chain)
                    try:
                        dictio_chainid[key_chain] = unique_chain_id
                    except:
                        print('Some error happened')
                        sys.exit(0)


    pdb_string, dict_conv = get_pdb_from_list_of_atoms(list_atoms, renumber=True, polyala=poliA, maintainCys=cys,
                                                       normalize=bfacnorm, sort_reference=True,
                                                       remove_non_res_hetatm=True, write_pdb=True,
                                                       dictio_chains=dictio_chainid,
                                                       path_output_pdb=model_file, predicted_model=predicted_model)

    # 3) Get the secondary and tertiary structure description with ALEPH
    graph_template, \
    oristru, \
    matrix, \
    cvs_list, \
    highd = ALEPH.annotate_pdb_model_with_aleph(pdb_model=model_file, weight="distance_avg", min_ah=min_alpha,
                                                min_bs=min_beta, write_pdb=True, strictness_ah=min_diff_ah,
                                                strictness_bs=min_diff_bs, peptide_length=peptide_length, is_model=False,
                                                only_reformat=True, compactness=predicted_model)

    graph_template = graph_template.vs.select(sstype_in=["ah","bs"]).subgraph()
    iterator_frag = ALEPH.get_all_fragments(graph_template)  # get only ah and bs, I will annotate coil later


    # 4) Open the pdb and populate the dictionary with the residues as keys and saving some information
    if len(Bio.PDB.Selection.unfold_entities(oristru, 'M')) > 1:
        print("\n Sorry, currently the use of NMR models is not supported in ARCIMBOLDO_SHREDDER")
        sys.exit(1)

    listres_oristru = Bio.PDB.Selection.unfold_entities(oristru, 'R')
    dict_oristru = {}
    for i, res in enumerate(listres_oristru):
        tuple_id = res.get_full_id()
        chain = tuple_id[2]
        nres = tuple_id[3][1]
        nameres = res.get_resname()

        dict_oristru[nres] = {'residue_object': res, 'ori_full_id': tuple_id, 'ori_chain': chain, 'ori_nres': nres,
                              'ori_nameres': nameres, 'ss_type_res': None, 'ss_reslist': None, 'ss_id_res': None,
                              'first_ref_cycle_group': None, 'second_ref_cycle_group': None,
                              'third_ref_cycle_group': None}

    # NOTE: At the moment, cycles correspond to:
    # first cycle: community clustering groups
    # second cycle: helices independent and beta strands as they were left by the community clustering
    # third cycle : same as second, but it is not performed by default

    # Get the total number of residues now, to compare later one how many are within secondary structure elements
    total_res = len(listres_oristru)
    print('\n Total number of residues in the input template is', total_res)

    # 5) Now we need to do the community clustering and populate the dict_oristru with the information
    if remove_coil and not gyre_preserve_chains:  # only then we want to annotate in groups
        vclust = ALEPH.get_community_clusters_one_step(algo=algorithm_community, graph_input_d=graph_template,
                                                       structure=oristru, pdb_search_in='', pathpdb='check',
                                                       write_pdb=True, pack_beta_sheet=pack_beta_community,
                                                       homogeneity=homogenity_community)
        if homogenity_community:
            vclust=vclust[0]
        results_community = ALEPH.get_dictionary_from_community_clusters(graph=graph_template, vclust=vclust,
                                                                         structure=oristru, writePDB=False,
                                                                         outputpath=None, header="",
                                                                         returnPDB=False, adding_t=0)

        list_groups = [results_community[key] for key in results_community.keys()]
        unique_community_groups = len(set(list_groups))

    print("\n Found ", len(iterator_frag), " secondary structure elements")
    listres_to_keep = []
    count_extra=0
    for i, dictio in enumerate(iterator_frag):
        # dictio is really a igraph.Vertex, but can be accessed as a dictionary. There is a dictio per fragment
        identifier_ss = dictio["sstype"] + str(i)
        if dictio["sstype"] == 'ah' or dictio["sstype"] == 'bs':
            filterlist = [tuple(x[:-1]) for x in dictio["reslist"]]
            print("\n     This is a ", dictio["sstype"], " of ", len(filterlist), " residues --> ",dictio['sequence'])
            if dictio["sstype"] == 'bs':
                print('Beta strands belong to beta sheet id ',dictio["sheet"])
            if (len(filterlist) < min_alpha and dictio["sstype"] == 'ah') or (len(filterlist) < min_beta and dictio["sstype"] == 'bs'):
                continue
            dictio["reslist"]=filterlist # we do not need sequence at this point, we can reduce it to its ids
            filterlist.sort()
            # Check if we need to elongate and how many residues
            if nres_extend != 0:  # only in the partial_coil case
                print('Entering the extension mode in partial_coil')
                # upwards
                first_res_id=filterlist[0][-1][1]
                elong_list_id=[]
                if first_res_id-nres_extend in dict_oristru:
                    #print '\nwe can extend upwards'
                    for indx in range(first_res_id-nres_extend,first_res_id):
                        # First check if we are entering in another secondary structure
                        if dict_oristru[indx]['ss_type_res']!=None:
                            #print 'We are not free to go, we are accessing another secondary structure'
                            break
                        dict_oristru[indx]['ss_type_res'] = dictio["sstype"]
                        dict_oristru[indx]['ss_id_res'] = identifier_ss
                        dict_oristru[indx]['first_ref_cycle_group']= dict_oristru[first_res_id]['first_ref_cycle_group']
                        dict_oristru[indx]["ori_full_id"]=('stru', dict_oristru[indx]["ori_full_id"][1],
                                                           dict_oristru[indx]["ori_full_id"][2],
                                                           dict_oristru[indx]["ori_full_id"][3])
                        filterlist.append(dict_oristru[indx]["ori_full_id"])
                        filterlist.sort()
                        elong_list_id.append(indx)
                    for indx in elong_list_id: # we want to iterate directly on the ids
                        dict_oristru[indx]['ss_reslist'] = filterlist
                        count_extra = count_extra + 1
                else:
                    continue

                # downwards
                last_res_id=filterlist[-1][-1][1]
                #print 'last_res_id',last_res_id
                if last_res_id+nres_extend in dict_oristru:
                    #print 'we can extend downwards'
                    #print 'filterlist before',filterlist
                    elong_list_id=[]
                    for indx in range(last_res_id+1,last_res_id+nres_extend+1):
                        #print 'indx',indx
                        #print 'dict_oristru[indx] before',dict_oristru[indx]
                        if dict_oristru[indx]['ss_type_res']!=None:
                            #print 'We are not free to go, we are accessing another secondary structure'
                            break
                        dict_oristru[indx]['ss_type_res'] = dictio["sstype"]
                        dict_oristru[indx]['ss_id_res'] = identifier_ss
                        dict_oristru[indx]['first_ref_cycle_group']= dict_oristru[last_res_id]['first_ref_cycle_group']
                        #print 'dict_oristru[indx] after',dict_oristru[indx]
                        #print 'dict_oristru[indx]["ori_full_id"] before',dict_oristru[indx]["ori_full_id"]
                        dict_oristru[indx]["ori_full_id"]=('stru', dict_oristru[indx]["ori_full_id"][1], dict_oristru[indx]["ori_full_id"][2], dict_oristru[indx]["ori_full_id"][3])
                        #print 'dict_oristru[indx]["ori_full_id"] after ',dict_oristru[indx]["ori_full_id"]
                        filterlist.append(dict_oristru[indx]["ori_full_id"])
                        filterlist.sort()
                        elong_list_id.append(indx)
                    #print 'filterlist after',filterlist
                    for indx in elong_list_id: # we want to iterate directly on the ids
                        dict_oristru[indx]['ss_reslist'] = filterlist
                        count_extra = count_extra + 1
                        #print 'Assigning new ss_reslist to ',indx
                        #print "dict_oristru[indx]['ss_reslist']",dict_oristru[indx]['ss_reslist']
                else:
                    continue
                    #print 'sorry, we can not extend downwards'

            listres_to_keep.extend(filterlist)

            if dictio["sstype"] == 'ah' or dictio["sstype"] == 'bs':
                for j in range(len(filterlist)):
                    key = filterlist[j][3][1]
                    dict_oristru[key]['ss_type_res'] = dictio["sstype"]
                    dict_oristru[key]['ss_id_res'] = identifier_ss
                    dict_oristru[key]['ss_reslist'] = filterlist


    # Annotate the dictionary by groups if community was performed
    if remove_coil and not gyre_preserve_chains:  # only then we want to annotate in groups
        # Check community results
        if unique_community_groups <= 0:
            print('\nWith current community clustering strategy there are no clusters returned, please check')
            sys.exit(1)

        # Get the correct ids to use and annotate groups for gyre and gimble
        list_used_ids=[]
        for key in sorted(results_community.keys()):
            tuple_id = key
            nres = tuple_id[3][1]
            dict_oristru[nres]['first_ref_cycle_group'] = results_community[key]
            ide = int(results_community[key].split("group")[1])
            if ide not in list_used_ids:
                list_used_ids.append(ide)
        list_used_ids.sort()
        last_used_id = list_used_ids[-1]

        # Perform the annotation for the other two levels
        number_id = last_used_id
        for i, dictio in enumerate(iterator_frag):
            print('\n Processing element ',i,' that corresponds to ',dictio)
            if dictio["sstype"] == 'ah':
                number_id += 1  # next group
                print('This is an alpha helix formed by ',dictio["reslist"])
                for j in range(len(dictio["reslist"])):
                    key = dictio["reslist"][j][3][1]
                    dict_oristru[key]['third_ref_cycle_group'] = 'group' + str(number_id)  # helices are treated independently in both second and third cycle grouping
                    dict_oristru[key]['second_ref_cycle_group'] = 'group' + str(number_id)
            elif dictio["sstype"] == 'bs':
                print('This is a beta strand formed by ',dictio["reslist"])
                number_id += 1  # next group
                for j in range(len(dictio["reslist"])):
                    key = dictio["reslist"][j][3][1]
                    dict_oristru[key]['third_ref_cycle_group'] = 'group' + str(number_id)  # one per beta strand
                    dict_oristru[key]['second_ref_cycle_group'] = dict_oristru[key]['first_ref_cycle_group'] # this will keep the same beta community as in the first grouping

    # 6) If remove_coil has been set, remove that residues from the model and the dictionary
    if remove_coil:
        for model in oristru.get_list():
            for chain in model.get_list():
                for residue in chain.get_list():
                    id_res = residue.get_full_id()
                    mini_id = residue.id
                    nres = mini_id[1]
                    if id_res not in listres_to_keep:
                        chain.detach_child(mini_id)
                        dict_oristru.pop(nres, None)
                    # NOTE CM: this was to avoid a problem if some residues had not been anotated.
                    # I think it is not needed, but in any case I should do it only if we are in the condition without
                    # user defined chains
                    if not gyre_preserve_chains:
                        if nres in dict_oristru and dict_oristru[nres]['first_ref_cycle_group'] is None:
                            chain.detach_child(mini_id)
                            dict_oristru.pop(nres, None)
    # NOTE CM: I am filtering out the coil residues from the graph, so I need this to annotate them
    # NOTE CM: At the moment they are annotated as single residues. Possibly doing it in continuous regions better
    else:
        count_coil = 0
        for model in oristru.get_list():
            for chain in model.get_list():
                for residue in chain.get_list():
                    id_res = residue.get_full_id()
                    mini_id = residue.id
                    nres = mini_id[1]
                    if dict_oristru[nres]['ss_type_res'] is None:
                        dict_oristru[nres]['ss_type_res'] = 'coil'
                        dict_oristru[nres]['ss_id_res'] = 'coil'+str(count_coil)
                        dict_oristru[nres]['ss_reslist'] = [id_res]
                        count_coil += 1


    if remove_coil and not gyre_preserve_chains:
        # Check if community clustering has assigned all residues in case we did it
        # Two possibilities: there is a problem with the distances or we extended the secondary structure elements
        if len(results_community) < len(dict_oristru.keys()) and nres_extend != 0:
            # it is normal that we don't have the same, check if the difference correspond to the number of added residues
            if len(results_community) + count_extra == len(dict_oristru.keys()):
                print('\n We are fine, we added extra ',count_extra,' residues ', ' Keep going ')
            else:
                print(colored("FATAL:", "red"), colored(
                "Community clustering did not asign a group to all residues, please consider modifying distances values for the community_clustering",
                'yellow'))
                sys.exit(1)
        elif len(results_community) < len(dict_oristru.keys()) and nres_extend == 0:
            if not gyre_preserve_chains: # in that case we are not going to consider community clustering groups
                print(colored("FATAL:", "red"), colored("Community clustering did not asign a group to all residues, please consider modifying distances values for the community_clustering",'yellow'))
                sys.exit(1)
        else:
            print('\nAll OK, Super!')

    # Check % of secondary structure
    ss_percentage=float(len(listres_to_keep))/total_res*100
    print('\n The percentage of secondary structure for this template is ',ss_percentage)
    if ss_percentage < 50:
        print(colored("""\nWARNING: With less than 50 per cent of secondary structure present in the template, it would be better to run ARCIMBOLDO_SHREDDER with the full template without removing the coil """,'yellow'))
    else:
        print(colored("\nMore than 50 per cent of the structure has secondary structure, continuing with the run", "magenta"))


    #Annotate domains
    dict_vclust = ALEPH.find_local_folds_in_the_graph(reference=model_file, strictness_ah=min_diff_ah,
                                                strictness_bs=min_diff_bs, peptide_length=peptide_length,
                                                write_pdb=False, min_ah=min_alpha, min_bs=min_beta, domains=True,
                                                compactness=predicted_model, write_graphml=False, write_graph_img=False)
    results_domain = {}
    for i, vclust in dict_vclust.items():
        results_domain = ALEPH.get_dictionary_from_community_clusters(graph=graph_template, vclust=vclust,
                                                            structure=oristru, writePDB=False,
                                                            outputpath=None, header="",
                                                            returnPDB=False, adding_t=0)
        for key in sorted(results_domain.keys()):
            nres = key[3][1]
            dict_oristru[nres]['domains_group_' + str(i)] = results_domain[key]


    # 5) Write the processed template to generate the models
    # Even the not gyre_preserve_chains condition has been checked before
    if not gyre_preserve_chains and remove_coil:
        print(colored("\nRemove coil has been set to on, and gyre and gimble will be performed according to automatic annotation", "magenta"))
        list_atoms = Bio.PDB.Selection.unfold_entities(oristru, 'A')
        list_atoms=sorted(list_atoms,key=lambda x:x.get_parent().get_full_id()[3][1]) # sort by residue number
        outpdbpath = os.path.join(current_directory, "shred_template.pdb")
        get_pdb_from_list_of_atoms(reference=list_atoms,renumber=False, uniqueChain=False, chainId='A',
                                   chainFragment=False, polyala=False, maintainCys=False,
                                   normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
                                   write_pdb=True, path_output_pdb=outpdbpath)
        model_file = os.path.abspath(os.path.join(current_directory, "shred_template.pdb"))

    if not gyre_preserve_chains and not remove_coil:
        print(colored("\nCoil has been left in the template model and no gyre or gimble refinement will be performed"
        , "magenta"))
        list_atoms = Bio.PDB.Selection.unfold_entities(oristru, 'A')
        list_atoms=sorted(list_atoms,key=lambda x:x.get_parent().get_full_id()[3][1]) # sort by residue number
        outpdbpath = os.path.join(current_directory, "shred_template.pdb")
        get_pdb_from_list_of_atoms(reference=list_atoms,renumber=False, uniqueChain=False, chainId='A',
                                   chainFragment=False, polyala=False, maintainCys=False,
                                   normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
                                   write_pdb=True, path_output_pdb=outpdbpath)
        model_file = os.path.abspath(os.path.join(current_directory, "shred_template.pdb"))

    if gyre_preserve_chains and remove_coil:
        print(colored("\nRemove coil has been set to on, and gyre and gimble will be performed according to user-given annotation", "magenta"))
        # save the model without the coil but with the chain it had at the beginning
        list_atoms = Bio.PDB.Selection.unfold_entities(oristru, 'A')
        outpdbpath = os.path.join(current_directory, "shred_template.pdb")
        list_atoms=sorted(list_atoms,key=lambda x:x.get_parent().get_full_id()[3][1]) # sort by residue number
        get_pdb_from_list_of_atoms(reference=list_atoms,renumber=False, uniqueChain=False, chainId="A",
                                   chainFragment=False, polyala=False, maintainCys=False,
                                   normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
                                   write_pdb=True, path_output_pdb=outpdbpath)
        model_file = os.path.abspath(os.path.join(current_directory, "shred_template.pdb"))
    if gyre_preserve_chains and not remove_coil:
        print(colored("\nThe coil has been left, and gyre and gimble will be performed according to user-given annotation", "magenta"))
        # save model with the coil but making sure non-annotated residues are not considered. Keep chains user-given
        list_atoms = Bio.PDB.Selection.unfold_entities(oristru, 'A')
        outpdbpath = os.path.join(current_directory, "shred_template.pdb")
        list_atoms = sorted(list_atoms, key=lambda x: x.get_parent().get_full_id()[3][1])  # sort by residue number
        get_pdb_from_list_of_atoms(reference=list_atoms,renumber=False, uniqueChain=False, chainId="A",
                                   chainFragment=False, polyala=False, maintainCys=False,
                                   normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
                                   write_pdb=True, path_output_pdb=outpdbpath)
        model_file = os.path.abspath(os.path.join(current_directory, "shred_template.pdb"))

    # Write the pdbs with their grouping levels to check
    if remove_coil and not gyre_preserve_chains:
        path_pdbfirst=os.path.join(current_directory,oristru_name+'_first_grouping_level.pdb')
        shutil.copy(model_file,path_pdbfirst)
        modify_chains_according_to_shredder_annotation(pdb=path_pdbfirst, dictio_annotation=dict_oristru,
                                       annotation_level='first_ref_cycle_group', output_path=current_directory)
        path_pdbsecond=os.path.join(current_directory,oristru_name+'_second_grouping_level.pdb')
        shutil.copy(model_file,path_pdbsecond)
        modify_chains_according_to_shredder_annotation(pdb=path_pdbsecond, dictio_annotation=dict_oristru,
                                       annotation_level='second_ref_cycle_group', output_path=current_directory)


    # 6) Get the distance matrix between the CA
    distance_matrix_CA, names_matrix = get_CA_distance_matrix(model_file)

    # 7) Save annotation in a pkl file
    save_annotation=open(os.path.join(current_directory,'annotated_template.pkl'),'wb')
    pickle.dump(dict_oristru,save_annotation,protocol=2)
    save_annotation.close()

    return model_file, dict_oristru, distance_matrix_CA, names_matrix, ss_percentage

def modify_dict_oristru(dict_oristru, model_file):
    model_name = os.path.basename(model_file)[:-4]
    structure_model = get_structure(model_name, model_file)
    listres_model = Bio.PDB.Selection.unfold_entities(structure_model, 'R')
    for res in listres_model:
        numres = res.get_full_id()[3][1]
        dict_oristru[numres]['residue_object'] = res

    return dict_oristru

def normalize_bfactors_of_pdb(pdbf, bf):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com
    :param pdbf:
    :type pdbf:
    :param bf:
    :type bf:
    :return:
    :rtype:
    """
    f = open(pdbf, "r")
    linee = f.readlines()
    f.close()

    f2 = open(pdbf, "w")
    bf = ("%.2f" % bf)
    if len(bf) <= 5:
        bf = ' ' + bf
    for linea in linee:
        if not linea.startswith("ATOM") and not linea.startswith("HETATM"):
            f2.write(linea)
        else:
            li = linea.split()
            lou = list(linea)
            lou[60] = bf[0]
            lou[61] = bf[1]
            lou[62] = bf[2]
            lou[63] = bf[3]
            lou[64] = bf[4]
            lou[65] = bf[5]
            mio = ''.join(lou)
            mio = mio.strip()
            f2.write(mio + "\n")
    f2.close()

def trim_sidechains_and_cysteines(pdb_model, poliA, cys):
    """ Modify a pdb to remove the sidechains and optionally keep cysteines

    :author: Claudia Millan
    :email: cmncri@ibmb.csic.es

    :param pdb_model:
    :type pdb_model:
    :param poliA:
    :type poliA:
    :param cys:
    :type cys:
    :return:
    :rtype:
    """
    # Prepare input pdb depending on choices and write it
    pdb_file = open(pdb_model, 'r')
    pdb_lines = pdb_file.readlines()
    pdb_file.close()
    pdb_medio = open(pdb_model, 'w')  # Overwrite the contents of the previous one
    for line in pdb_lines:
        if not line.startswith("ANISO") and not line.startswith("ATOM") and not line.startswith("HETATM"):
            pdb_medio.write(line)
        elif line.startswith("ATOM") or line.startswith("HETATM"):
            parts = line.split()
            type_res = parts[3]
            if type_res.endswith('HOH'):
                continue  # go to next line
            if line.startswith("HETATM"):
                if type_res not in ['MSE','SEP','TPO','MIR']: # NOTE: should I include more non-standard residues?
                    continue # go to next line
            if poliA == False and cys == False:
                pdb_medio.write(line)  # We do not need to do any selection, just write
            elif poliA == True and cys == True:  # Then we have to save only poliA but mantain the cysteines
                list_items = list(line)
                atom = ''.join(list_items[13:16])
                type_at = atom.strip()
                if type_at in ["CA", "CB", "N", "C", "O", "SG"]:
                    pdb_medio.write(line)
            elif poliA == True and cys == False:  # Plain poliala
                list_items = list(line)
                atom = ''.join(list_items[13:16])
                type_at = atom.strip()
                if type_at in ["CA", "CB", "N", "C", "O"]:
                    pdb_medio.write(line)
            elif poliA == False and cys == True:  # This is redundant, if it has sidechains it will have its cysteins, Same as first option
                pdb_medio.write(line)
    pdb_medio.close()

def modify_chains_according_to_shredder_annotation(pdb, dictio_annotation, annotation_level, output_path):
    """ Given a pdb or a list of pdbs, using the annotation dictionary, uses one of the annotation levels
    and produce and rewrites the pdbs with that chain definition

    :author: Claudia Millan
    :email: cmncri@ibmb.csic.es

    :param pdb: paths to the pdbs to modify
    :type pdb: str or list of str
    :param dictio_annotation: following format one key per each residue inside
           419: {'residue_object': <Residue ARG het=  resseq=419 icode= >, 'ss_type_res': 'bs', 'ori_nameres': 'ARG',
           'ori_nres': 419, 'ss_reslist': [('1hdh_0_0', 0, 'A', (' ', 419, ' ')), ('1hdh_0_0', 0, 'A', (' ', 420, ' ')),
           ('1hdh_0_0', 0, 'A', (' ', 421, ' ')), ('1hdh_0_0', 0, 'A', (' ', 420, ' ')), ('1hdh_0_0', 0, 'A', (' ', 421, ' '))]
           , 'first_ref_cycle_group': 'group0', 'third_ref_cycle_group': 'group21', 'second_ref_cycle_group': 'group0',
           'ss_id_res': 'bs36', 'ori_chain': 'A', 'ori_full_id': ('1hdh_0_0', 0, 'A', (' ', 419, ' '))}
    :type dictio_annotation: dict
    :param annotation_level: can be: 'third_ref_cycle_group','second_ref_cycle_group','first_ref_cycle_group'
    :type annotation_level: str
    :param output_path: path where the pdb(s) with the changes in the annotation must be written
    :type output_path: str
    :return:
    :rtype:
    """
    dictio_chainid = {}
    if not isinstance(pdb, list):
        pdb = [pdb]
    for i, pdb_file in enumerate(pdb):
        structure = get_structure(name=os.path.basename(pdb_file[:-4]), pdbf=pdb_file)
        for model in structure:
            for chain in model:
                for residue in chain:
                    key_chain = residue.get_full_id()  # that is the key for the dictio_chains dictionary
                    key_annotation = key_chain[3][1]
                    group = int(dictio_annotation[key_annotation][annotation_level][5:])
                    indx_group = group
                    try:
                        # the chain id must be different per each different group!
                        dictio_chainid[key_chain] = list_id[indx_group]
                    except:
                        print('There are too many groups defined, there are not any more possible ids')
                        sys.exit(0)
        outputpdb_path = os.path.join(output_path, os.path.basename(pdb_file))
        pdb_file_atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')
        pdb_file_atoms = sorted(pdb_file_atoms, key=lambda x:x.get_parent().get_full_id()[3][1])  # sort by res number
        get_pdb_from_list_of_atoms(reference=pdb_file_atoms, path_output_pdb=outputpdb_path, dictio_chains=dictio_chainid,
                                   normalize=False, sort_reference=True, renumber=False, uniqueChain=False,
                                   polyala=False, maintainCys=False, write_pdb=True)

def get_CA_distance_matrix(pdb_model):
    """

    :param pdb_model:
    :type pdb_model:
    :return:
    :rtype:
    """
    structure = get_structure(name=pdb_model[:-4], pdbf=pdb_model)
    list_for_matrix=[]
    names_matrix_dict={}
    list_residues = Bio.PDB.Selection.unfold_entities(structure,'R')
    list_CA_atoms = [ residue['CA'] for _,residue in enumerate(list_residues)]
    for i in range(len(list_CA_atoms)):
        list_for_matrix.append([])
        for j in range(len(list_CA_atoms)):
            #print 'i,list_CA_atoms[i].get_full_id()[3][1],j,list_CA_atoms[j].get_full_id()[3][1]',i,list_CA_atoms[i].get_full_id()[3][1],j,list_CA_atoms[j].get_full_id()[3][1]
            id_first=list_CA_atoms[i].get_full_id()[3][1]
            id_second=list_CA_atoms[j].get_full_id()[3][1]
            list_for_matrix[i].append(list_CA_atoms[i]-list_CA_atoms[j])
            if id_first in names_matrix_dict:
                names_matrix_dict[id_first][id_second]=[i,j]
            else:
                names_matrix_dict[id_first]={id_second: [i,j]}
    distance_CA_matrix=numpy.array(list_for_matrix)
    numpy.set_printoptions(precision=3)
    #print '\n',distance_CA_matrix
    #print '\n', names_matrix_dict
    return distance_CA_matrix,names_matrix_dict

def elongate_extremities(pdb_model, dictio_template, list_distances_full, res_to_complete, target_size, min_ah=7, min_bs=4):
    """ Extends the models generated by shredder_spheres by looking at the extremities of the fragments.

    :param pdb_model:
    :type pdb_model: str
    :param dictio_template:
    :type dictio_template:
    :param list_distances_full:
    :type list_distances_full:
    :param res_to_complete:
    :type res_to_complete:
    :param min_ah:
    :type min_ah:
    :param min_bs:
    :type min_bs:
    :return:
    :rtype:
    """

    trials_limit = 0
    list_already_elongated = []
    structure = get_structure(name=pdb_model[:-4], pdbf=pdb_model)
    list_residues = sorted(Bio.PDB.Selection.unfold_entities(structure, 'R'),
                           key=lambda x: x.get_full_id()[3][1:])  # list of residues objects sorted by id
    list_initial_nres = [resi.get_full_id()[3][1] for resi in list_residues]
    list_initial_idss = list(set([dictio_template[nres]['ss_id_res'] for nres in list_initial_nres]))

    while trials_limit < 201:  # number of times you will perform the iterative process as a maximum

        trials_limit = trials_limit + 1

        # Reread the PDB
        structure = get_structure(name=pdb_model[:-4], pdbf=pdb_model)
        list_resi_sorted = sorted(Bio.PDB.Selection.unfold_entities(structure, 'R'), key=lambda x: x.get_full_id()[3][1:])
        list_distances = []
        list_nres_elongations = []
        list_removal = []

        # Get the extremities of the continues stretches on the model
        list_extremities = []
        continuous = []
        continuous_fraglist = []
        for i, resi in enumerate(list_resi_sorted):
            if i < len(list_resi_sorted) - 1:
                check = check_continuity(res1=resi, res2=list_resi_sorted[i + 1])
                if not check:
                    list_extremities.append(('end', resi))
                    list_extremities.append(('start', list_resi_sorted[i + 1]))
                    if resi.get_full_id()[3][1] not in continuous:
                        continuous.append(resi.get_full_id()[3][1])
                    continuous_fraglist.append(continuous)
                    continuous = []
                    if i == len(list_resi_sorted) - 2:
                        continuous_fraglist.append([list_resi_sorted[i + 1].get_full_id()[3][1]])
                else:
                    if resi.get_full_id()[3][1] not in continuous:
                        continuous.append(resi.get_full_id()[3][1])
                    if list_resi_sorted[i + 1].get_full_id()[3][1] not in continuous:
                        continuous.append(list_resi_sorted[i + 1].get_full_id()[3][1])
            else:
                break
        if continuous != [] and (continuous not in continuous_fraglist):
            continuous_fraglist.append(continuous)

        # You need to add the two extremes that are so by definition, that is, the first and the last residue
        list_extremities.append(('end', list_resi_sorted[-1]))
        list_extremities.append(('start', list_resi_sorted[0]))

        continuous_fraglist.sort(key=lambda x: len(x))

        list_idsizetuple = []
        for stretch in continuous_fraglist:
            list_residues_ss = []
            for residue in stretch:
                ss_ident = dictio_template[residue]['ss_id_res']
                list_residues_ss.append((residue, ss_ident))

            list_groups = [list(group) for key, group in itertools.groupby(list_residues_ss, operator.itemgetter(1))]

            for current in list_groups:
                current_group_ss = current[0][1]  # the first one gives us the key
                list_resi = [ele[0] for ele in current]
                list_idsizetuple.append((current_group_ss, list_resi))

        list_idsizetuple.sort(key=lambda x: len(x[1]))

        if res_to_complete == 0:
            if len(list_idsizetuple[0][1]) >= min_ah and list_idsizetuple[0][0].startswith('ah'):
                break
            if len(list_idsizetuple[0][1]) >= min_bs and list_idsizetuple[0][0].startswith('bs'):
                break
            if len(list_idsizetuple[0][1]) >= 1 and list_idsizetuple[0][0].startswith('coil'):
                break
            for i in range(len(list_idsizetuple)):
                cont_section_ss = list_idsizetuple[i]
                ss_ident_section = dictio_template[cont_section_ss[1][0]]['ss_id_res']
                list_res = [x[-1][1] for x in dictio_template[cont_section_ss[1][0]]['ss_reslist']]
                if ss_ident_section.startswith('ah') and len(cont_section_ss[1]) < min_ah:
                    remain = min_ah - len(cont_section_ss[1])
                elif ss_ident_section.startswith('bs') and len(cont_section_ss[1]) < min_bs:
                    remain = min_bs - len(cont_section_ss[1])
                # rest is in common
                not_yet = [val for val in list_res if val not in cont_section_ss[1]]
                for l in range(1, len(list_idsizetuple)):
                    longest_ss_id = list_idsizetuple[-l][0]
                    longest_ss_res = list_idsizetuple[-l][1]
                    if longest_ss_id == ss_ident_section:
                        continue
                    else:
                        break
                # check that longest_ss does not contain already residues we want to add
                if len(set(not_yet).intersection(set(longest_ss_res))) != 0:  # THIS SHOULD NOT HAPPEN!
                    print('There was a problem, please report to bugs-borges@ibmb.csic.es')
                    sys.exit(0)
                else:
                    list_removal.extend(longest_ss_res[-remain:])
                for index, _ in enumerate(not_yet):
                    if index == len(list_removal):
                        break
                    list_nres_elongations.append(not_yet[index])
                break
        elif res_to_complete < 0:
            print('Something very wrong happened, please report to bugs-borges@ibmb.csic.es')
            sys.exit(1)
        else:
            if trials_limit == 201:
                if res_to_complete > 10:
                    # Then I want to exclude this model
                    return dictio_template, False, list_nres_elongations
                else:
                    # Then I am fine, just return OK
                    return dictio_template, True, list_nres_elongations
        # Retrieve the distance from the CA defined in center to all the extremities
        for _, residue_tuple in enumerate(list_extremities):
            residue = residue_tuple[1]
            tag_ext = residue_tuple[0]
            nres = residue.get_full_id()[3][1]
            for i, _ in enumerate(list_distances_full):
                if list_distances_full[i][0] == nres:
                    list_distances.append((list_distances_full[i][1], residue.get_full_id(), tag_ext))
        list_distances.sort(key=lambda x: x[0])

        # Check that the extremes can be elongated
        for i, _ in enumerate(list_distances):
            tag_ext = list_distances[i][2]
            resi_full_id = list_distances[i][1]
            nres = (resi_full_id[3][1])
            # modify reslist to check the numerical id only
            filter_reslist = [x[-1] for x in dictio_template[nres]['ss_reslist']]
            preext_ind = (filter_reslist).index(resi_full_id[-1])
            # two possibilities for extension, depending on n or cterminal
            if preext_ind == 0 and tag_ext == 'start':  # we can't extend below
                continue
            elif preext_ind == len(filter_reslist) - 1 and tag_ext == 'end':  # we can't extend after
                continue
            else:
                if tag_ext == 'end':
                    sup_add = filter_reslist[preext_ind + 1][1]
                    if (sup_add not in list_already_elongated) and (sup_add not in list_initial_nres):
                        list_nres_elongations.append(sup_add)
                elif tag_ext == 'start':
                    inf_add = filter_reslist[preext_ind - 1][1]
                    if (inf_add not in list_already_elongated) and (inf_add not in list_initial_nres):
                        list_nres_elongations.append(inf_add)

        # If there are no points for elongation we need to go to another secondary structure element
        if len(list_nres_elongations) == 0:
            if res_to_complete < min(min_ah,min_bs): # If adding a new one is going to make it incomplete
                break
            # Otherwhise search for the next ss element that I can add
            for i, _ in enumerate(list_distances_full):
                if (list_distances_full[i][0] not in list_already_elongated) and \
                        (list_distances_full[i][0] not in list_initial_nres) and \
                        (dictio_template[list_distances_full[i][0]]['ss_id_res'] not in list_initial_idss):
                    next_res = list_distances_full[i][0]
                    list_nres_elongations.extend([next_res])
                    break
                else:
                    continue

        # Remove what needs to be removed
        if len(list_removal) > 0:
            for model in structure.get_list():
                for chain in model.get_list():
                    for residue in chain.get_list():
                        nres_res = residue.get_full_id()[3][1]
                        if nres_res in list_removal:
                            try:
                                index_res = list_already_elongated.index(nres_res)
                                del (list_already_elongated[index_res])
                            except:
                                pass
                            try:
                                index_res2 = list_initial_nres.index(nres_res)
                                del (list_initial_nres[index_res2])
                            except:
                                pass
                            mini_id = residue.id
                            chain.detach_child(mini_id)
                            res_to_complete = res_to_complete + 1
            # flush list_removal
            list_removal = []

        # Perform the actual elongation of the model
        for nres in list_nres_elongations:
            addresobj = dictio_template[nres]['residue_object']
            id_chain_where_to_add = dictio_template[nres]['ori_chain']
            if (res_to_complete > 0):  # to be sure we don't add more than we need
                try:
                    structure[0][id_chain_where_to_add].add(addresobj)
                except Exception as e:
                    if str(e).endswith('defined twice'):  # Then the residue is already in the structure
                        continue
                    elif str(e).replace("'","") == str(id_chain_where_to_add):
                        structure[0].add(Bio.PDB.Chain.Chain(id_chain_where_to_add))
                        structure[0][id_chain_where_to_add].add(addresobj)
                        pass
                    else:
                        print('Something wrong happened with the pdb construction: ', e)
                        sys.exit(1)
                res_to_complete = res_to_complete - 1
                list_already_elongated.append(nres)


        # Save the model before the next iteration
        list_atoms_sorted = sorted(Bio.PDB.Selection.unfold_entities(structure, 'A'),
                                   key=lambda x: x.get_parent().get_full_id()[3][1:])  # sort by residue number
        #         get_pdb_from_list_of_atoms(reference=new_list_atoms, renumber=False, polyala=False, maintainCys=False,
        #                           normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
        #                           write_pdb=True, path_output_pdb=pdbmodel_path)
        get_pdb_from_list_of_atoms(reference=list_atoms_sorted, path_output_pdb=pdb_model, write_pdb=True,
                                   renumber=False, uniqueChain=False, polyala=False, maintainCys=False,
                                   normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={})
        # Flush the list of the elongations
        list_nres_elongations = []

    return dictio_template, True, list_nres_elongations

def filter_models_by_coordinates(path_folder):
    """ Check a folder of pdbs, if they are identical in coordinates, remove the copies and leave one representative.

    :param path_folder:
    :type path_folder:
    :return:
    :rtype:
    """
    # NOTE CM: this could be improved a lot by using something faster and safer
    redundant = []
    complete_list = os.listdir(path_folder)
    dict_coord = {}
    for __index, name in enumerate(complete_list):
        structure = get_structure(name=name[:-4], pdbf=os.path.join(path_folder, name))
        dict_coord[name] = sorted([list(a.get_coord()) for a in Bio.PDB.Selection.unfold_entities(structure, 'A')],
                                  key=lambda x: (x[0], x[1], x[2]))
    # Now we can start to compare
    for i, name1 in enumerate(complete_list):
        coord1 = dict_coord[name1]
        for j, name2 in enumerate(complete_list):
            if j <= i:
                continue
            coord2 = dict_coord[name2]
            if len(coord1) != len(coord2):
                continue
            if all([coord1[k] == coord2[k] for k in range(len(coord1))]):
                redundant.append(name2)
    # print "len(redundant)",len(set(redundant))
    for __index, redundant_model in enumerate(set(redundant)):
        os.remove(os.path.join(path_folder, redundant_model))

def shredder_spheres(working_directory, namedir, pdb_model, dictio_template, target_size, dist_matrix, convNamesMatrix,
                     min_ah=7, min_bs=4, step=1, list_centers=[]):
    """ Generates a set of compact models.

    :param working_directory: current working directory in the SHREDDER run
    :type working_directory: str
    :param namedir: output folder
    :type namedir: str
    :param pdb_model: template to extract the models from
    :type pdb_model: str
    :param dictio_template: annotation of the pdb given in pdb_model in terms of secondary and tertiary structure
    :type dictio_template: dict
    :param target_size: size that the models should have
    :type target_size: int
    :param dist_matrix:
    :type dist_matrix:
    :param convNamesMatrix:
    :type convNamesMatrix:
    :param min_ah:
    :type min_ah:
    :param min_bs:
    :type min_bs:
    :param step: step to traverse the sequence
    :type step: int
    :param list_centers: [(name_frag(str),center(list of floats with x,y,z coord of atom)),...]
                        if not empty, only the CA in the list will be used for model generation
    :type list_centers: list of tuples
    :return:
    :rtype:
    """

    # Recognize the path and create the library folder
    current_wd = working_directory
    wd_library = os.path.join(current_wd, namedir)
    if not os.path.exists(wd_library):
        os.mkdir(wd_library)
    elif os.path.exists(wd_library):
        shutil.rmtree(wd_library)
        os.mkdir(wd_library)

    # Generate input for the search
    print('######################################################################################')
    print('GENERATING MODELS FROM STARTING TEMPLATE AT ', pdb_model)
    print('######################################################################################')

    # Get the list with all CAs in the template
    structure = get_structure(name=pdb_model[:-4], pdbf=pdb_model)
    list_all_CA_atoms = [residue['CA']
                         for residue in Bio.PDB.Selection.unfold_entities(structure, 'R')
                         if residue.has_id('CA')]


    # Check the list_centers option and act accordingly
    if not list_centers:  # Then we want to use ALL the centers defined by the number of residues and step
        for i in range(0, len(list_all_CA_atoms), step):
            nfragtag = str((list_all_CA_atoms[i].get_full_id())[3][1])
            center = list_all_CA_atoms[i].get_coord()
            list_centers.append((nfragtag, center))

    log_distance_models = open('distance_models_shredder.txt', 'w')
    for i, _ in enumerate(list_centers):

        # Identify which is the residue from we will be generating the model
        name_frag = list_centers[i][0]
        print('\n \n \n \n **=====******======******======****======********===********===***********===******')
        print('\n Processing model centered in ', name_frag)

        # Obtain the list of pairwise distances between this residue and all the rest
        sorted_dist_list = get_sorted_distance_list_to_CA(dist_matrix, convNamesMatrix, name_frag)

        # Prepare the lists in which we will perform the model generation
        residues_on_list = [[sorted_dist_list[x][0], 'off'] for x, _ in enumerate(sorted_dist_list)]
        new_list_sort = [[sorted_dist_list[x][0],
                          sorted_dist_list[x][1],
                          dictio_template[sorted_dist_list[x][0]]['ss_id_res']]
                         for x, _ in enumerate(sorted_dist_list)]
        only_nres = [new_list_sort[ni][0] for ni, _ in enumerate(new_list_sort)]
        only_nres_sorted = sorted(only_nres)


        # Start the model generation
        residues_on = 0

        for ires, resi in enumerate(new_list_sort):

            if residues_on > target_size:
                print('\n We passed the target size, reducing the model')
                added.sort(key=operator.itemgetter(1))
                resremove = abs(target_size - residues_on)
                count = 0
                for ind, element in enumerate(added):
                    if count == resremove:
                        break
                    indx = only_nres.index(element[1])
                    residues_on_list[indx][1] = 'off'
                    residues_on -= 1
                    count += 1

            if residues_on == target_size:
                print('\n We reached exactly the target size')
                break

            if target_size - max(min_ah,min_bs) < residues_on <= target_size:
                print('\n We almost reached the target size')
                break

            if residues_on_list[ires][1] == 'on':  # we skip it, it was already selected
                continue
            else:
                resi_ss = resi[2]
                prev_resi_id = resi[0]
                if resi_ss.startswith('bs'):  # limit to minimal size of beta strands
                    limit_ss = min_bs
                elif resi_ss.startswith('ah'):  # limit to minimal size of alpha helices
                    limit_ss = min_ah
                elif resi_ss.startswith('coil'):  # limit to minimal size of coil
                    limit_ss = 1
                # Check what other things we can add
                added = []
                for inext in range(ires, len(new_list_sort)):
                    if len(added) >= limit_ss:
                        break
                    current_resi_id = new_list_sort[inext][0]
                    current_resi_ss = new_list_sort[inext][2]
                    if current_resi_ss == resi_ss:
                        # If the closest residue of the same ss is not on, set it
                        if residues_on_list[inext][1] == 'off':
                            added.append((inext, current_resi_id))
                            residues_on_list[inext][1] = 'on'
                            residues_on = residues_on + 1
                            # Add the remaining residues in between these two
                            current_ind = only_nres_sorted.index(current_resi_id)
                            prev_ind = only_nres_sorted.index(prev_resi_id)
                            if abs(current_ind - prev_ind) > 1:
                                mini = min(current_ind, prev_ind)
                                maxi = max(current_ind, prev_ind)
                                list_cutind = only_nres_sorted[mini + 1:maxi]
                                for ele in list_cutind:
                                    iord = only_nres.index(ele)
                                    if residues_on_list[iord][1] == 'off':
                                        residues_on_list[iord][1] = 'on'
                                        residues_on = residues_on + 1
                                        added.append((iord, ele))

                        prev_resi_id = current_resi_id


        # Write the pdb with the residues that are set to on
        to_save_nres = [ele[0] for ele in residues_on_list if ele[1] == 'on']
        structure = get_structure(name=pdb_model[:-4], pdbf=pdb_model)
        io = Bio.PDB.PDBIO()
        io.set_structure(structure)
        for model in structure.get_list():
            for chain in model.get_list():
                for residue in chain.get_list():
                    residue_nres = residue.get_full_id()[3][1]
                    if not (residue_nres in to_save_nres):
                        chain.detach_child(residue.id)
        pdbmodel_path = wd_library + "/" + 'frag' + name_frag + '_0_0.pdb'
        new_list_atoms = Bio.PDB.Selection.unfold_entities(structure, 'A')
        new_list_atoms = sorted(new_list_atoms,
                                key=lambda x: x.get_parent().get_full_id()[3][1])  # sort by residue number

        #write_pdb_file_from_list_of_atoms(new_list_atoms, pdbmodel_path, dictio_chains={}, renumber=False, uniqueChain=False)
        get_pdb_from_list_of_atoms(reference=new_list_atoms, renumber=False, polyala=False, maintainCys=False,
                                   normalize=False, sort_reference=True, remove_non_res_hetatm=False, dictio_chains={},
                                   write_pdb=True, path_output_pdb=pdbmodel_path)

        # Check the new current model size
        list_nres_elongations=[]
        list_ca = [ atom for atom in new_list_atoms if atom.id == 'CA']
        size = len(list_ca)
        print('\n     Size before elongation', size)
        diff_res = target_size - size
        if diff_res > 0:
            # Note CM: including a condition to avoid doing elongation if we are already really close to target
            if diff_res > min(min_ah,min_bs):  # we can test different values here
                print('\n     We are going to extend the secondary structure elements of the model by ',diff_res)
                dictio_template, elong_bool, list_nres_elongations  = elongate_extremities(pdb_model=pdbmodel_path, dictio_template=dictio_template,
                                                                   list_distances_full=sorted_dist_list, target_size=target_size,
                                                                   res_to_complete=diff_res, min_ah=min_ah, min_bs=min_bs)
                if elong_bool == True:  # elongation went OK, we want to keep this model
                    print('\n This model was correctly elongated, we will keep it')
                else:
                    print('============================================================================')
                    print('\n Something went wrong with the elongation of this model, we will skip it')
                    print('============================================================================')
                    os.remove(pdbmodel_path)
            else:
                print('\n     This model is really close to size required, skipping elongation')

        else:
            print('\n     This model has already the size required')


        if len(list_nres_elongations)>0:
            distances_list = [ele for ele in sorted_dist_list if (ele[0] in to_save_nres) or (ele[0] in list_nres_elongations)]
            distances_only = [ele[1] for ele in distances_list]
        else:
            distances_list = [ ele for ele in sorted_dist_list if ele[0] in to_save_nres ]
        distances_only=[ele[1] for ele in distances_list]
        log_distance_models.write(pdbmodel_path+'\t'+str(max(distances_only))+'\n')

    del log_distance_models

    # Filter the redundancy that might remain anyway
    # TODO CM: Filter using a different method
    print("\n********************************************************************************")
    print('Before filtering there are ', len(os.listdir(wd_library)), ' models')
    filter_models_by_coordinates(wd_library)
    print('After filtering there are ', len(os.listdir(wd_library)), ' models')
    print("********************************************************************************\n")

    return dictio_template


def get_sorted_distance_list_to_CA(dist_matrix,convNamesMatrix,center_name):
    """

    :param dist_matrix:
    :type dist_matrix:
    :param convNamesMatrix:
    :type convNamesMatrix:
    :param center_name:
    :type center_name:
    :return:
    :rtype:
    """
    #print ' I need to get all distances to residue ',center_name
    #print 'check where to search in the matrix, convNamesMatrix[int(center_name)]',convNamesMatrix[int(center_name)]
    list_distances=[]
    key_center=int(center_name)
    list_possible_keys=convNamesMatrix.keys()
    for _,key in enumerate(list_possible_keys):
        index_i=convNamesMatrix[key_center][key][0]
        index_j=convNamesMatrix[key_center][key][1]
        # Get the values with the correct index
        value_distance=dist_matrix[index_i][index_j]
        list_distances.append((key,value_distance))
    list_distances.sort(key = lambda x: x[1])
    #return list_distances[1:]
    return list_distances

def write_csv_node_file(g, reference,pdbsearchin, min_diff_ah, min_diff_bs):

    o = open("write_csv_node_file.log", 'w')

    ss_record = pdbsearchin.split('\n')
    ss_record_dict = {}

    def catch_substring(a, b, line):
        try:
            return line[a:b].strip()
        except:
            return ""

    for i, element in enumerate(ss_record):
        if ss_record[i] != '':
            ss_record[i]=(element.split())
            o.write(str(ss_record[i]))

            if element[0:6].strip() == 'HELIX':
                ss_record_dict[catch_substring(19, 20, element) + '_' + catch_substring(21, 25, element) + '_' + catch_substring(33, 37, element)] =\
                    {'ss': element[0:6].strip(), 'ss_position': '1', 'chain': catch_substring(19, 20, element),
                                       'res_start': catch_substring(21, 25, element), 'res_end': catch_substring(33, 37, element), 'ss_direction': '0', 'ss_id': catch_substring(7, 10, element)}
                ss_id = int(catch_substring(7, 10, element))

            elif element[0:6].strip() == 'SHEET':
                if int(catch_substring(7, 10, element)) == 1:
                    ss_id += 1
                ss_record_dict[catch_substring(21, 22, element) + '_' + catch_substring(22, 26, element) + '_' +
                               catch_substring(33, 37, element)] ={'ss': element[0:6].strip(), 'ss_position': catch_substring(7, 10, element), 'chain': catch_substring(21, 22, element),
                                       'res_start': catch_substring(22, 26, element), 'res_end': catch_substring(33, 37, element), 'ss_direction': catch_substring(38, 40, element), 'ss_id':ss_id}

            else:
                o.write("\nThis is not a conventional Helix/Sheet record")
                o.write(str(element[0:6].strip()))
        else:
            ss_record.remove(ss_record[i])

    #o.write(str(ss_record_dict))

    o.write("\n===========================================================================================================")

    f = open('csv_node_file.csv', 'w')
    f.write("identifier,pdb_id,model,chain,res_start,res_end,ah_strict,bs_strict,ss_type,ss_id,ss_position,ss_direction,sequence")
    for frag in g.vs:
        #print(frag)
        identifier = (os.path.basename(reference)[:-4]+'_' + str(frag['reslist'][0][1])+'_'+
                      str(frag['reslist'][0][2]) + '_' + str(frag['reslist'][0][3][1])+'_'+str(frag['reslist'][-1][3][1])
                      +'_'+ str(min_diff_ah) +'_'+ str(min_diff_bs))
        if frag['sstype'] == 'coil': #Si el fragmento del grafo es coil
            pass

        elif str(frag['reslist'][0][2]) + '_' + str(frag['reslist'][0][3][1])+'_'+str(frag['reslist'][-1][3][1]) not in ss_record_dict: #Hay algçun tipo de error en comparar el grafo y el diccionario
            o.write("\nThere is a dismatch in your identifiers\n")
            o.write(str(frag['reslist'][0][2]) + '_' + str(frag['reslist'][0][3][1]) + '_' + str(frag['reslist'][-1][3][1]))
            #print(frag)
            print((str(frag['reslist'][0][2])) + '_' + str(frag['reslist'][0][3][1])+'_'+str(frag['reslist'][-1][3][1]))

        elif str(frag['reslist'][0][2]) + '_' + str(frag['reslist'][0][3][1])+'_'+str(frag['reslist'][-1][3][1]) in ss_record_dict: #En el caso que se pueda comparar el fragmetno delgrafo y el diccionario
            key = str(frag['reslist'][0][2]) + '_' + str(frag['reslist'][0][3][1])+'_'+str(frag['reslist'][-1][3][1])
            #print(key)
            f.write('\n')
            f.write(identifier + ',' + os.path.basename(reference)[:-4] + ',' + str(frag['reslist'][0][1]) + ',' + str(frag['reslist'][0][2]) + ',' + str(frag['reslist'][0][3][1]) + ','
                    + str(frag['reslist'][-1][3][1]) + ',' + str(min_diff_ah) + ','+ str(min_diff_bs) + ',' + str(frag['sstype']) + ',' + str(ss_record_dict[key]['ss_id']) + ',' +
                    str(ss_record_dict[key]['ss_position']) + ',' + str(ss_record_dict[key]['ss_direction']) + ',' + str(frag['sequence']))

            o.write((identifier + ',' + os.path.basename(reference)[:-4] + ',' + str(frag['reslist'][0][1]) + ',' +
                    str(frag['reslist'][0][2]) + ',' + str(frag['reslist'][0][3][1]) + ',' +
                    str(frag['reslist'][-1][3][1]) + ',' + str(min_diff_ah) + ','+ str(min_diff_bs) +',' + str(frag['sstype']) + ',' + str(ss_record_dict[key]['ss_id']) + ',' +
                    str(ss_record_dict[key]['ss_position']) + ',' + str(ss_record_dict[key]['ss_direction'])
                    + ',' + str(frag['sequence'])))
            o.write('\n')

            ss_record_dict[frag.index]=str(identifier) #Añadir una llave con la correspondencia entre el indice del nodo y el identificador

        else: #Error no esperado
            o.write("Not expected condition")
    f.close()
    o.close()
    return(ss_record_dict)

def write_csv_edge_file (g, ss_record_dict):
    #STILL IN DEVELOPING
    f = open('csv_edge_file.csv', 'w')
    f.write("node1_id,node2_id,weight\n")

    o = open("write_csv_node_file.log", 'w')

    for frag in g.es:
        try:
            f.write(str(ss_record_dict[frag.source]) + ',' + str(ss_record_dict[frag.target]) + ',' + str(frag['weight']) + '\n')
        except:
            o.write("Somethhing happened during the edge processing")
            try:
                o.write(frag)
            except:
                print("No frag")

    f.close()
    o.close()

def get_atoms_list(pdb):
    """ Generate a list with the atoms from a pdb file.

    :author: Ana del Rocío Medina Bernal from Bioinformatics.py
    :email: ambcri@ibmb.csic.es
    :param pdb: path to the pdb file
    :type pdb: str
    :return list_atoms:
    :rtype list_atoms: list
    """
    f = open(pdb, "r")
    lines = f.readlines()
    f.close()
    list_atoms = []

    for line in lines:
        lispli = line.split()
        if len(lispli) > 0 and lispli[0] in ["ATOM", "HETATM"]:
            list_atoms.append(lispli)

    return list_atoms

def change_chain(pdb, chain, atom_list=["ATOM  ", "ANISOU", "HETATM", "TER   "]):
    """ Generate a pdb file with the Chain modified

    :author: Ana del Rocío Medina Bernal from Bioinformatics.py
    :email: ambcri@ibmb.csic.es
    :param pdb: path to the pdb file
    :type pdb: str
    :param chain: chain name
    :type chain: str
    :param atom_list: pdb lines to change
    :type atom_list: str
    :return out: pdb file modified
    :rtype list_atoms: str
    """

    f = open(pdb, "r")
    lines = f.readlines()
    f.close()
    out = ""
    for line in lines:
        # only look at records indicated by atom_list
        if line[0:6] not in atom_list:
            if not line.startswith("END"):
                out += line + "\n"
            continue
        # Grab only residues belonging to chain
        out += line[:21] + chain + line[22:] + "\n"
    return out

def pdb_offset(pdb, offset):
    """
    Adds an offset to the residue column of a pdb file without touching anything
    else.

    :author: Ana del Rocío Medina Bernal from Bioinformatics.py
    :email: ambcri@ibmb.csic.es
    :param pdb: path to the pdb file
    :type pdb: str
    :param offset: chain name
    :type offset: int
    :return out: pdb file modified
    :rtype list_atoms: str
    """

    # Read in the pdb file
    f = open(pdb,'r')
    lines = f.readlines()
    f.close()

    out = []
    for line in lines:
        # For and ATOM record, update residue number
        if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
            num = offset + int(line[22:26])
            out.append("%s%4i%s" % (line[0:22], num, line[26:]))
        else:
            out.append(line)

    return "\n".join(out)

def get_number_of_residues(pdb):
    """
       Counts the number of residues of a pdb file

       :author: Ana del Rocío Medina Bernal from Bioinformatics.py
       :email: ambcri@ibmb.csic.es
       :param pdb: path to the pdb file
       :type pdb: str
       :return number: number of residues
       :rtype list_atoms: int
       """
    stru = get_structure("ref", pdb)
    number = 0
    for model in stru.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                if residue.has_id("CA") and residue.has_id("C") and residue.has_id("O") and residue.has_id("N"):
                    number += 1
    return number

def get_ideal_helices_from_lenghts(list_length, pdb, reversed=False):
    helix_list = []
    helipdb = io.StringIO(str(pdb))
    stru = get_structure("helix", helipdb)
    for helix_length in list_length:
        if isinstance(helix_length, list) or isinstance(helix_length, tuple):
            helix_list.append(get_ideal_helices_from_lenghts(helix_length, pdb, reversed=reversed))
        else:
            cont = 0
            lips = []
            allcont = 0
            for model in stru:
                for chain in model:
                    for res in chain:
                        if reversed and allcont < 70-helix_length:
                            allcont += 1
                            continue

                        allcont += 1
                        if cont < helix_length:
                            lips += res.get_list()
                            cont += 1
                        else:
                            break
            helipdb = get_pdb_from_list_of_atoms(lips)
            helix_list.append(helipdb)
    return helix_list

#FROM Bioinformatics.py

#TODO: Rename this function to complain the format standard: generate_pdb
# used in SELSLIB2
def generatePDB(pdbf, resi_list, filepath, trim_to_polyala=True):
    root, basename = os.path.split(filepath)
    if not os.path.exists(root):
        os.makedirs(root)
    structure = get_structure("ref", pdbf)
    atoms_list = []
    for model in structure:
        for chain in model:
            for resi in chain:
                if resi.get_id()[1] in resi_list:
                    if trim_to_polyala:
                        atoms_list += [resi["CA"], resi["C"], resi["O"], resi["N"]]
                        if resi.has_id("CB"):
                            atoms_list.append(resi["CB"])
                    else:
                        for atom in resi:
                            atoms_list.append(atom)

    pdball = get_pdb_from_list_of_atoms(atoms_list, chainFragment=True)[0]
    f = open(filepath, "w")
    f.write(pdball)
    f.close()

#TODO: Rename this function to complain the format standard: generate_pdb_omitting
# used in SELSLIB2
def generatePDBomitting(pdbf, resi_list, filepath, trim_to_polyala=True):
    """

    :param pdbf:
    :type pdbf:
    :param resi_list:
    :type resi_list:
    :param filepath:
    :type filepath:
    :param trim_to_polyala:
    :type trim_to_polyala:
    :return:
    :rtype:
    """
    root, basename = os.path.split(filepath)
    if not os.path.exists(root):
        os.makedirs(root)
    structure = get_structure("ref", pdbf)
    atoms_list = []
    for model in structure:
        for chain in model:
            for resi in chain:
                if resi.get_id()[1] not in resi_list:
                    if trim_to_polyala:
                        atoms_list += [resi["CA"], resi["C"], resi["O"], resi["N"]]
                        if resi.has_id("CB"):
                            atoms_list.append(resi["CB"])
                    else:
                        for atom in resi:
                            atoms_list.append(atom)
    try:
        pdball = get_pdb_from_list_of_atoms(atoms_list, chainFragment=True)[0]
        f = open(filepath, "w")
        f.write(pdball)
        f.close()
    except:
        print("Cannot generate", filepath, "from")
        print("resi_list", resi_list)

#TODO: Rename this function to complain the format standard: sequential_renumber_list_of_pdbs
# used in SELSLIB2
def sequentialRenumberListOfPdbs(pdblist,keep_chains=False):
    out = []
    chain = "A"
    for pdb in pdblist:
        for line in pdb.splitlines():
            # For and ATOM record, update residue number
            if line[0:6] == "ATOM  " or line[0:6] == "TER   ":
                if keep_chains:
                    chain=line[21]
                num = 0 + int(line[22:26])
                epse = "%s%4i%s" % (line[0:21] + chain, num, line[26:])
                out.append(epse)
            else:
                out.append(line)
        chain = "" + str(chr((ord(chain) + 1)))
    return "\n".join(out)

#TODO: SELSLIB2 is using this function although I think we should move on from this
@SystemUtility.deprecated("This function is not yhet supported and should be not used anymore")
def __STDCentroidKParameter(data, v_cycles, treshJump, criteria="mean", minKappa=None, maxKappa=None, oneByone=False):
    narr = numpy.array(data)
    whitened = scipy.cluster.vq.whiten(narr)
    # whitened = vq.whiten(narr)
    valori = []

    kappa = 0
    lastkappa = 0
    start = (numpy.sqrt(len(whitened) / 2) / 2)

    if minKappa == None:
        kappa = start
    else:
        kappa = minKappa

    if maxKappa == None:
        lastkappa = start * 2 * 2
    else:
        lastkappa = maxKappa

    lastkappa = len(whitened)

    i = int(kappa)
    if i <= 0:
        i = 1

    centroids = []
    while True:
        if i > lastkappa:
            break
        print("trying kappa", i)
        sys.stdout.flush()
        avg_crossv = 0.0
        avg_sqd = 0.0
        kappa_step = 0
        allcent = []
        alldist = []
        for v in range(v_cycles):
            centroids, distorsion = scipy.cluster.vq.kmeans(whitened, i, iter=20)
            allcent.append(centroids)
            alldist.append(distorsion)

        alldist = numpy.array(alldist)
        minInd = numpy.argmin(alldist)
        distorsion = None
        if criteria == "max":
            distorsion = numpy.max(alldist)
        elif criteria == "min":
            distorsion = numpy.min(alldist)
        elif criteria == "mean":
            distorsion = numpy.mean(alldist)
        else:
            print("Error!!!, criteria, ", criteria, " not recognized!")
            raise Exception("Criteria function argument not recognized!")

        centroids = allcent[minInd]
        # print "----",centroids
        # print distorsion
        # groups, labels = vq.kmeans2(training,i,iter=20,minit="points")
        print("done...")
        # D_k = [scipy.spatial.distance.cdist(whitened, [cent], 'euclidean') for cent in groups]
        # cIdx = [numpy.argmin(D,axis=1) for D in D_k]
        # dist = [numpy.min(D,axis=1) for D in D_k]
        # print "D_k",D_k
        # print ";;;",dist
        # avgWithinSS = [sum(d)/whitened.shape[0] for d in dist]
        # print ",,,,,",avgWithinSS

        """
        sum_sqd = 0
        for centroid in groups:
            print "......",centroid
            print labels
            print "-----",numpy.std(centroid)
            sum_sqd += numpy.std(centroid)
        avg_sqd = sum_sqd / len(groups)
        """
        valori.append([distorsion, i])
        print(i, distorsion)
        kappa = i
        if len(valori) > 1:
            jump = (valori[-2])[0] - (valori[-1])[0]
            print("difference", (valori[-2])[0] - (valori[-1])[0])
            if not oneByone and distorsion >= 3:
                i += 150
                print("Next Kappa will be", i, "increased by", 150)
            elif not oneByone and not oneByone and distorsion >= 2:
                i += 100
                print("Next Kappa will be", i, "increased by", 100)
            elif not oneByone and distorsion >= 1:
                i += 80
                print("Next Kappa will be", i, "increased by", 80)
            elif not oneByone and distorsion >= 0.7:
                i += 70
                print("Next Kappa will be", i, "increased by", 70)
            elif not oneByone and distorsion >= 0.6:
                i += 60
                print("Next Kappa will be", i, "increased by", 60)
            elif not oneByone and distorsion >= 0.5:
                i += 50
                print("Next Kappa will be", i, "increased by", 50)
            elif not oneByone and distorsion >= 0.4:
                i += 40
                print("Next Kappa will be", i, "increased by", 40)
            elif not oneByone and distorsion >= 0.3:
                i += 30
                print("Next Kappa will be", i, "increased by", 30)
            elif not oneByone and distorsion >= 0.2:
                i += 20
                print("Next Kappa will be", i, "increased by", 20)
            elif not oneByone and distorsion >= 0.1:
                i += 10
                print("Next Kappa will be", i, "increased by", 10)
            else:
                i += 1
                print("Next Kappa will be", i, "increased by", 1)

            if jump <= -0.6:
                # kappa -= 1
                break

            if len(valori) > 2 and abs(jump) <= abs(treshJump):
                # kappa -= 1
                break
        else:
            i += 1
            print()

    code, distance = scipy.cluster.vq.vq(whitened, centroids)
    # print "***",code
    return kappa, code

#TODO: Used in ARCIMBOLDO_BORGES and SELSLIB2
# Remove the call to getSuperimp and instead call a specific function in ALEPH to superpose big structures

def getListCA(name, pdbf, mode, algorithm="biopython", DicParameters=None, backbone=False, listmodel=None,
              allInList=False, minResInChain=None, superpose_exclude=1, pos=[], clusterCYS=False):
    """

    :param name:
    :type name:
    :param pdbf:
    :type pdbf:
    :param mode: Can be 'DB','PDB','PRECOMPUTED','PDBSTRINGBM', 'PDBSTRING', 'PDBSTRINGBM_RESIDUES_CONSERVED'
    :type mode:
    :param algorithm:
    :type algorithm:
    :param DicParameters:
    :type DicParameters:
    :param backbone:
    :type backbone:
    :param listmodel:
    :type listmodel:
    :param allInList:
    :type allInList:
    :param minResInChain:
    :type minResInChain:
    :param superpose_exclude:
    :type superpose_exclude:
    :param pos:
    :type pos:
    :param clusterCYS:
    :type clusterCYS:
    :return:
    :rtype:
    """
    if mode == "DB" and DicParameters != None:
        corresponding = (pdbf.split("/"))[-1]
        pdbid, model, idSolution = corresponding.split("_")
        idSolution, ext = idSolution.split(".")
        nameExp = DicParameters["nameExperiment"]
        (DicParameters, db) = SystemUtility.requestConnectionDatabase(DicParameters)
        cur = db.cursor()
        cur.execute(
            "SELECT backbone,pda FROM " + nameExp + "_SOLUTIONS WHERE pdbid = '" + pdbid + "' AND model = " + model + " AND IdSolution = " + idSolution)
        nodo = cur.fetchone()
        back, pda = nodo
        reference = pickle.loads(back)
        cur.close()
        db.close()
        return reference
    elif mode == "PDB" and listmodel == None:
        structure = get_structure(name, pdbf)
        reference = []
        lineallCA = []
        lineallwCB = []
        for model in structure.get_list():
            for chain in model.get_list():
                if minResInChain != None and len(chain.get_list()) < minResInChain:
                    continue
                for residue in chain.get_list():
                    if residue.has_id("CA"):
                        atom = residue["CA"]
                        reference.append(atom)
                        if allInList:
                            if residue.has_id("CB"):
                                lineallCA.append([residue["CA"], residue["C"], residue["O"], residue["N"], residue["CB"]])
                                lineallwCB.append([residue["CA"], residue["C"], residue["O"], residue["N"]])
                            else:
                                lineallwCB.append([residue["CA"], residue["C"], residue["O"], residue["N"]])
                                lineallCA.append([residue["CA"], residue["C"], residue["O"], residue["N"]])

                    if backbone:
                        # print pdbf,residue.get_full_id()
                        if residue.has_id("CB"):
                            reference.append(residue["CB"])
                        if residue.has_id("C"):
                            reference.append(residue["C"])
                        if residue.has_id("O"):
                            reference.append(residue["O"])
                        if residue.has_id("N"):
                            reference.append(residue["N"])
        withoutcb = []
        for at in reference:
            if at.get_id() != "CB":
                withoutcb.append(at)
        if allInList:
            return [[lineallCA, lineallwCB]]
        else:
            return [[reference, withoutcb]]
    else:
        allfrags = []
        if mode == "PDB":
            allfrags = getFragmentsListWithoutCVS(pdbf)
        elif mode == "PRECOMPUTED":
            allfrags = pdbf
        elif mode.startswith("PDBSTRINGBM"):
            try:
                stru_template = get_structure("stru",io.StringIO(str(pdbf)))
                list_resi = Bio.PDB.Selection.unfold_entities(stru_template, 'R')
                list_resi_sorted=sorted(list_resi,key=lambda x:x.get_full_id()[3][1:])
                allfrags = [list_resi_sorted]
            except:
                print(sys.exc_info())
                traceback.print_exc(file=sys.stdout)
        elif mode == "PDBSTRING":
            # allfrags = getFragmentsListWithoutCVS(pdbf,notafile=True)
            # tuplone = getFragmentListFromPDBUsingAllAtoms(cStringIO.StringIO(pdbf),False)
            strucd = get_structure("ddds", io.StringIO(str(pdbf)))
            # asde = range(len(listmodel))
            apt = []
            for model in strucd.get_list():
                for chain in model.get_list():
                    trtrtr = sorted(chain.get_unpacked_list(), key=lambda x: x.get_id()[1])

                    # DONE: Change to Continuity by C-N distance
                    for r in range(0, len(trtrtr) - 1):
                        resi1 = trtrtr[r]
                        resi2 = trtrtr[r + 1]
                        # print resi1.get_id()
                        # print resi2.get_id()
                        if checkContinuity(resi1, resi2):
                            # if resi2.get_id()[1] != resi1.get_id()[1]+1:
                            #   print "Really contigous??",resi2.get_id()[1],resi1.get_id()[1]
                            apt.append(resi1)
                            # print "r",r,"r+1",r+1,"len(trtrtr)-2",len(trtrtr)-2,len(trtrtr)
                            if r == len(trtrtr) - 2:
                                apt.append(resi2)
                        else:
                            apt.append(resi1)
                            allfrags.append(apt)
                            apt = []
                            if r == len(trtrtr) - 2:
                                apt.append(resi2)

                    if len(apt) > 0:
                        allfrags.append(apt)
                        apt = []
                    """
                    for residue in sorted(chain.get_unpacked_list(),__resBioOrder2):
                            #print residue
                            apt.append(residue)
                            if len(apt) == listmodel[asde[0]]["fragLength"]:
                                    allfrags.append(deepcopy(apt))
                                    apt = []
                                    del asde[0]
                    """

        #print "??????????????????????????????????????",len(allfrags)
        if mode == "PDBSTRINGBM":
            listCombi = itertools.permutations(allfrags)
        elif mode == "PDBSTRINGBM_RESIDUES_CONSERVED":
            listCombi = [allfrags]
        else:
            listCombi = [allfrags]

        listValid = []
        if listmodel != None:
            """
            for fra in listmodel:
                 print "frag",fra["fragLength"]
            print "============================="
            for combi in listCombi:
                for crazy in combi:
                        print "crazy",len(crazy)
                print "-----------------"

            listCombi = itertools.permutations(allfrags)
            """
            lio = []
            for combi in listCombi:
                lio.append(combi)
            listCombi = lio

            for combi in listCombi:
                valid = False
                for ui in range(len(listmodel)):
                    frag1 = listmodel[ui]
                    frag2 = combi[ui]
                    # print "Comparing",len(frag2), frag1["fragLength"]
                    if len(frag2) != frag1["fragLength"]:
                        valid = False
                        break
                    else:
                        valid = True
                # print "is valid",valid

                if valid:
                    torange = 1
                    if algorithm in ["biopython-core", "nigels-core2"]:
                        torange = superpose_exclude
                    fragli = [[] for _ in range(len(combi))]
                    for qw in range(torange):
                        for pids in range(len(combi)):
                            it = combi[pids]
                            itans = __sublisttrim(it, qw)
                            listof = []
                            for ita in itans:
                                reference = []
                                for rt in range(len(ita)):
                                    residue = ita[rt]
                                    # if residue.has_id("CA"):
                                    # if (not clusterCYS or (clusterCYS and residue.get_resname() in ["CYS","cys","Cys"] and rt in pos)) and residue.has_id("CA"):
                                    if not clusterCYS and residue.has_id("CA"):
                                        atom = residue["CA"]
                                        reference.append(atom)
                                    # if backbone:
                                    # if (not clusterCYS or (clusterCYS and residue.get_resname() in ["CYS","cys","Cys"] and rt in pos)) and backbone:
                                    if not clusterCYS and backbone:
                                        # print pdbf,residue.get_full_id()
                                        if residue.has_id("CB"):
                                            reference.append(residue["CB"])
                                        if residue.has_id("C"):
                                            reference.append(residue["C"])
                                        if residue.has_id("O"):
                                            reference.append(residue["O"])
                                        if residue.has_id("N"):
                                            reference.append(residue["N"])
                                    if clusterCYS and residue.get_resname() in ["CYS", "cys", "Cys"] and rt in pos:
                                        if residue.has_id("CB"):
                                            reference.append(residue["CB"])
                                        if residue.has_id("SG"):
                                            reference.append(residue["SG"])

                                listof.append((reference, len(ita)))
                            fragli[pids] += listof

                    for com in itertools.product(*fragli):
                        reference = []
                        lenfrags = []
                        for fragitem in com:
                            frag, lenfra = fragitem
                            reference += frag
                            lenfrags.append(lenfra)

                        withoutcb = []
                        for at in reference:
                            if at.get_id() != "CB":
                                withoutcb.append(at)
                        listValid.append([reference, withoutcb, lenfrags])
        else:
            for combi in listCombi:
                torange = 1
                if algorithm in ["biopython-core", "nigels-core2"]:
                    torange = superpose_exclude
                fragli = [[] for _ in range(len(combi))]
                for qw in range(torange):
                    for pids in range(len(combi)):
                        it = combi[pids]
                        itans = __sublisttrim(it, qw)
                        listof = []
                        for ita in itans:
                            reference = []
                            for rt in range(len(ita)):
                                residue = ita[rt]
                                # if residue.has_id("CA"):
                                # if (not clusterCYS or (clusterCYS and residue.get_resname() in ["CYS","cys","Cys"] and rt in pos)) and residue.has_id("CA"):
                                if not clusterCYS and residue.has_id("CA"):
                                    atom = residue["CA"]
                                    reference.append(atom)
                                # if backbone:
                                # if (not clusterCYS or (clusterCYS and residue.get_resname() in ["CYS","cys","Cys"] and rt in pos)) and backbone:
                                if not clusterCYS and backbone:
                                    if residue.has_id("CB"):
                                        reference.append(residue["CB"])
                                    if residue.has_id("C"):
                                        reference.append(residue["C"])
                                    if residue.has_id("O"):
                                        reference.append(residue["O"])
                                    if residue.has_id("N"):
                                        reference.append(residue["N"])
                                if clusterCYS and residue.get_resname() in ["CYS", "cys", "Cys"] and rt in pos:
                                    if residue.has_id("CB"):
                                        reference.append(residue["CB"])
                                    if residue.has_id("SG"):
                                        reference.append(residue["SG"])
                            listof.append((reference, len(ita)))
                        fragli[pids] += listof

                #print 'len product fragli',len(list(product(*fragli)))
                #print "analysis",[map(lambda x: x[1], a) for a in fragli]

                for com in itertools.product(*fragli):
                    reference = []
                    lenfrags = []
                    for fragitem in com:
                        frag, lenfra = fragitem
                        reference += frag
                        lenfrags.append(lenfra)
                    withoutcb = []
                    for at in reference:
                        if at.get_id() != "CB":
                            withoutcb.append(at)
                    listValid.append([reference, withoutcb, lenfrags])

                    # print "Validi are",len(listValid)
                    # if len(listValid) == 0:
                    #               for fra in listmodel:
                    #                 print "frag_refe",fra["fragLength"]
                    #               for fra in allfrags:
                    #                print "frag_solu",len(fra)
        return listValid
        # print "You need to use PDB and DB mode with or without listmodel, but if it is PRECOMPUTED  mode you need listmodel set"

#TODO: This function is needed by getSuperimp but must be deleted with it
def __sublisttrim(fragment, delres):
    allist = []
    listfra1 = fragment[delres:]
    listfra2 = fragment[:-1 * (delres)]
    # print "Sublist trim function:",len(fragment),"delres is",delres
    # print listfra1
    # print "===============**************=================="
    # print listfra2
    # print "listafrag1: ",len(listfra1)
    # print "listafrag2: ",len(listfra2)

    if len(listfra1) > 0:
        allist.append(listfra1)
    if len(listfra2) > 0:
        allist.append(listfra2)

    for i in range(1, delres):
        resto = -1 * (delres - i)
        listfra3 = fragment[i:resto]
        # print "listafrag_i: ",len(listfra3)
        if len(listfra3) > 0:
            allist.append(listfra3)
    return allist

#TODO: This function is needed by getSuperimp but must be deleted with it
def getListAllAtoms(name, pdbf):
    structure = None
    try:
        structure = get_structure(name, pdbf)
    except:
        structure = get_structure(name, io.StringIO(str(pdbf)))

    reference = []
    for model in structure.get_list():
        for chain in model.get_list():
            for residue in chain.get_list():
                for atom in residue:
                    reference.append(atom)
    return structure, reference

#TODO: This function is needed by getSuperimp but must be deleted with it
def writePDBFromListOfAtom(reference, outputFilename, dictio_chains={}, renumber=False, uniqueChain=False):
    '''Write a pdb from a list of atoms.
    Input:
    - reference: string with input path
    - outputFilename: string with output path
    - dictio_chains: dictionary, keys are residues full id, and values are the new chain id
    - renumber: boolean
    - uniqueChain: boolean '''
    pdbString = ""
    numero = 0
    previous = None
    numea = 1
    for item in reference:
        orig_atom_num = item.get_serial_number()
        hetfield, resseq, icode = item.get_parent().get_id()
        segid = item.get_parent().get_segid()
        if renumber:
            hetfield = " "
            orig_atom_num = numea
            numea += 1
            if previous == None or resseq != previous:
                numero += 1
                previous = resseq
                resseq = numero
            elif previous != None:
                resseq = numero
            icode = " "
        if uniqueChain:
            chain_id = "A"
        if item.get_parent().get_full_id() in dictio_chains:
            chain_id = dictio_chains[item.get_parent().get_full_id()]
        else:
            chain_id = item.get_parent().get_parent().get_id()
        resname = item.get_parent().get_resname()
        element = item.get_name()
        pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id)
    flu = open(outputFilename, "w")
    flu.write(pdbString)
    flu.close()

#TODO: This function is needed by getSuperimp but must be deleted with it
def checkContinuity(res1, res2):
    resaN = None
    prevResC = None
    try:
        resaN = res2["N"]
        prevResC = res1["C"]
    except:
        return False

    # print resaC.get_coord()
    # print prevResN.get_coord()
    resaNX = float(resaN.get_coord()[0])
    resaNY = float(resaN.get_coord()[1])
    resaNZ = float(resaN.get_coord()[2])
    prevResCX = float(prevResC.get_coord()[0])
    prevResCY = float(prevResC.get_coord()[1])
    prevResCZ = float(prevResC.get_coord()[2])
    checkCont = numpy.sqrt(((resaNX - prevResCX) ** 2) + ((resaNY - prevResCY) ** 2) + ((resaNZ - prevResCZ) ** 2))
    # print resan
    # print prev_model,prev_chain,prevRes
    # print checkCont
    if checkCont <= 1.5:
        return True
    else:
        return False

@SystemUtility.deprecated('Use new get_pdb_from_list_of_atoms')
def getPDBFromListOfAtom(reference, renumber=False, uniqueChain=False, chainId="A", chainFragment=False, diffchain=None,
                         polyala=True, maintainCys=False):
    pdbString = ""
    numero = 1
    resn = {}
    nur = 1
    lastRes = None
    prevChain = 0
    lich = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "L", "M", "N", "O", "P", "Q", "R", "S", "T", "U", "V", "Z",
            "J", "K", "X", "Y", "W", "a", "b", "c", "d", "e", "f", "g", "h", "i", "l", "m", "n", "o", "p", "q", "r",
            "s", "t", "u", "v", "z", "j", "k", "x", "y", "w", "1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
    # if chainFragment or uniqueChain:
    #   renumber = True
    erel = []
    if not polyala:
        for item in reference:
            res = item.get_parent()
            for atk in res:
                erel.append(atk)
        reference = erel
    elif maintainCys:
        for item in reference:
            res = item.get_parent()
            if res.get_resname() in ["CYS", "cys", "Cys"]:
                for atk in res:
                    erel.append(atk)
            else:
                erel.append(item)
        reference = erel

    for item in reference:
        orig_atom_num = item.get_serial_number()
        hetfield, resseq, icode = item.get_parent().get_id()
        segid = item.get_parent().get_segid()
        chain_id = item.get_parent().get_parent().get_id()
        hetfield = " "
        if (resseq, chain_id) not in resn.keys():
            if lastRes != None:
                # print "Checking Continuity",lastRes.get_parent().get_full_id(),item.get_parent().get_full_id(),checkContinuity(lastRes.get_parent(),item.get_parent())
                # print "Checking Continuity",item.get_parent().get_full_id(),lastRes.get_parent().get_full_id(),checkContinuity(item.get_parent(),lastRes.get_parent())
                # print lich[prevChain]
                # print
                if not checkContinuity(lastRes.get_parent(), item.get_parent()) and not checkContinuity(
                        item.get_parent(), lastRes.get_parent()):
                    if renumber:
                        nur += 10
                    if chainFragment:
                        prevChain += 1
            new_chain_id = chain_id
            if uniqueChain:
                new_chain_id = chainId
            elif chainFragment:
                new_chain_id = lich[prevChain]
            if renumber:
                resn[(resseq, chain_id)] = (nur, new_chain_id)
            else:
                resn[(resseq, chain_id)] = (resseq, new_chain_id)
            lastRes = item
            nur += 1
        tuplr = resn[(resseq, chain_id)]
        resseq = tuplr[0]
        chain_id = tuplr[1]
        icode = " "
        orig_atom_num = numero
        numero += 1

        resname = item.get_parent().get_resname()
        element = item.get_name()
        pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id)

    if diffchain != None and len(diffchain) > 0:
        prevChain += 1
        lastRes = None
        for item in diffchain:
            orig_atom_num = item.get_serial_number()
            hetfield, resseq, icode = item.get_parent().get_id()
            segid = item.get_parent().get_segid()
            chain_id = item.get_parent().get_parent().get_id()
            hetfield = " "
            if (resseq, chain_id) not in resn.keys():
                if lastRes != None:
                    if not checkContinuity(lastRes.get_parent(), item.get_parent()) and not checkContinuity(
                            item.get_parent(), lastRes.get_parent()):
                        if renumber:
                            nur += 10
                        if chainFragment:
                            prevChain += 1
                new_chain_id = chain_id
                if uniqueChain:
                    new_chain_id = chainId
                elif chainFragment:
                    new_chain_id = lich[prevChain]

                if renumber:
                    resn[(resseq, chain_id)] = (nur, new_chain_id)
                else:
                    resn[(resseq, chain_id)] = (resseq, new_chain_id)
                lastRes = item
                nur += 1
            tuplr = resn[(resseq, chain_id)]
            resseq = tuplr[0]
            chain_id = tuplr[1]
            icode = " "
            orig_atom_num = numero
            numero += 1

            resname = item.get_parent().get_resname()
            element = item.get_name()
            pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id)

    return pdbString, resn

#TODO: getRMSDfromCommonCore is used in SELSLIB2 and should be renamed into get_rmsd_from_common_core
def getRMSDfromCommonCore(dicc_1, dicc_2):
    stru1_refe = dicc_1["reference"]
    stru1_cmp = dicc_1["compare"]
    stru2_refe = dicc_2["reference"]
    stru2_cmp = dicc_2["compare"]

    core = {}
    list1_CA = []
    list2_CA = []
    listc_CA = []

    for key in dicc_1:
        if key not in ["reference", "compare"]:
            if key in dicc_2:
                core[key] = (dicc_1[key][0], dicc_2[key][0], dicc_1[key][1], dicc_1[key][2] - dicc_2[key][2])
                # print "---",dicc_1[key][0]
                list1_CA.append(dicc_1[key][0]["CA"].get_coord())
                list2_CA.append(dicc_2[key][0]["CA"].get_coord())
                listc_CA.append(dicc_1[key][1]["CA"].get_coord())

    if len(core.keys()) > 0:
        ref1 = numpy.array(list1_CA)
        ref2 = numpy.array(list2_CA)
        comp = numpy.array(listc_CA)

        transf, rmsd_list = csb.bio.utils.fit_wellordered(ref1, comp, n_iter=None, full_output=True)
        if len(rmsd_list) > 0:
            rmsd_1 = rmsd_list[-1][1]
        else:
            rmsd_1 = -1
        transf, rmsd_list = csb.bio.utils.fit_wellordered(ref2, comp, n_iter=None, full_output=True)
        if len(rmsd_list) > 0:
            rmsd_2 = rmsd_list[-1][1]
        else:
            rmsd_2 = -2
    else:
        # for key in dicc_1:
        #       print "+++",key
        # for key in dicc_2:
        #       print "---",key
        # quit()
        rmsd_1 = -1
        rmsd_2 = -1

    return (rmsd_1, rmsd_2, core)

#TODO: this is used in SELSLIB2, but in a function executePicasso that is not called anymore... maybe it is the case to remove this
def getRMSD(referenceFile, compareFile, mode, DicParameters=None, algorithm="biopython", backbone=False, listmodel=None,
            doNotMove=False):
    if isinstance(referenceFile, str):
        listValiRef = getListCA("referencia", referenceFile, mode, DicParameters=DicParameters, backbone=backbone,
                                listmodel=listmodel)
        reference = (listValiRef[0])[1]
    elif isinstance(referenceFile, list):
        reference = referenceFile
        listValiRef = [[reference, reference]]

    if isinstance(compareFile, str):
        listValiCmp = getListCA("compare", compareFile, mode, DicParameters=DicParameters, backbone=backbone,
                                listmodel=listmodel)
        compare = (listValiCmp[0])[1]
    elif isinstance(compareFile, list):
        compare = compareFile
        listValiCmp = [[compare, compare]]

    nref = len(reference)
    ncom = len(compare)
    # print "nref",nref,"ncom",ncom
    # print "algorithm",algorithm
    pdball = ""

    if algorithm != "biopython":
        nameficheref = "./ref_" + str(datetime.datetime.now()) + ".pdb"
        writePDBFromListOfAtom(reference, nameficheref)

        namefichecmp = "./cmp_" + str(datetime.datetime.now()) + ".pdb"
        writePDBFromListOfAtom(compare, namefichecmp)

    # for atom in reference:
    #       print atom.get_full_id()
    # print "===="
    # for atom in compare:
    #       print atom.get_full_id()

    try:
        if algorithm == "biopython":
            # print mode, compareFile
            if not isinstance(compareFile, str) and mode != "PDBSTRING" and isinstance(compareFile[0], list):
                listValiCmp = getListValidFromTwoList(reference, compare)

            best_rmsd = 1000000000000
            best_ncom = -1
            for tren in range(len(listValiCmp)):
                no_cmp, compare, lfr = listValiCmp[tren]
                ncom = len(compare)
                # print "nref",nref,"ncom",ncom
                rmsd = None
                if doNotMove:
                    r = []
                    c = []
                    for atm in reference:
                        r.append(atm.get_coord())
                    for atm in compare:
                        c.append(atm.get_coord())
                    refe = numpy.array(r)
                    compa = numpy.array(c)
                    print(len(reference), len(compare))
                    diff = refe - compa
                    l = refe.shape[0]
                    rmsd = numpy.sqrt(sum(sum(diff * diff)) / l)
                else:
                    super_imposer = Bio.PDB.Superimposer()
                    super_imposer.set_atoms(reference, compare)
                    rmsd = super_imposer.rms

                if rmsd <= best_rmsd:
                    best_rmsd = rmsd
                    best_ncom = ncom

            rmsd = best_rmsd
            ncom = best_ncom
            # print "--",rmsd,ncom
        elif algorithm == "pymol_align":
            namefichepml = "./pml_" + str(datetime.datetime.now()) + ".pml"
            flu = open(namefichepml, "w")
            flu.write("load " + nameficheref + ", \"ref\"\n")
            flu.write("load " + namefichecmp + ", \"cmp\"\n")
            flu.write("align cmp and name ca, ref and name ca\n")
            flu.close()
            p = subprocess.Popen(['pymol', "-c", namefichepml], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            listout = out.split("\n")
            rmsd = -100
            for line in listout:
                listline = line.split()
                if line.startswith(" Executive: RMS ="):
                    rmsd = float(listline[3])
            os.remove(namefichepml)
        elif algorithm == "pymol_cealign":
            namefichepml = "./pml_" + str(datetime.datetime.now()) + ".pml"
            flu = open(namefichepml, "w")
            flu.write("load " + nameficheref + ", \"ref\"\n")
            flu.write("load " + namefichecmp + ", \"cmp\"\n")
            flu.write("cealign ref and name ca, cmp and name ca\n")
            flu.close()
            p = subprocess.Popen(['pymol', "-c", namefichepml], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            listout = out.split("\n")
            rmsd = -100
            for line in listout:
                listline = line.split()
                if line.startswith("RMSD"):
                    rmsd = float(listline[1])
            os.remove(namefichepml)
        elif algorithm == "pymol_super":
            namefichepml = "./pml_" + str(datetime.datetime.now()) + ".pml"
            flu = open(namefichepml, "w")
            flu.write("load " + nameficheref + ", \"ref\"\n")
            flu.write("load " + namefichecmp + ", \"cmp\"\n")
            flu.write("super cmp and name ca, ref and name ca\n")
            flu.close()
            p = subprocess.Popen(['pymol', "-c", namefichepml], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            listout = out.split("\n")
            rmsd = -100
            for line in listout:
                listline = line.split()
                if line.startswith(" Executive: RMS ="):
                    rmsd = float(listline[3])
            os.remove(namefichepml)
        if algorithm != "biopython":
            os.remove(nameficheref)
            os.remove(namefichecmp)
    except:
        rmsd = -100
        print(sys.exc_info())
        traceback.print_exc(file=sys.stdout)

    return (rmsd, nref, ncom)

@SystemUtility.deprecated("This function is not anymore supported and will be substituted soon")
def getSuperimp(referenceFile, compareFile, mode, DicParameters=None, algorithm="biopython", backbone=False,
                getDictioCorresp=False, superpose_exclude=1, listmodel=None, allAtomsList=None, allAtomsModel=None,
                n_iter=None, pos=[], clusterCYS=False, minmaxrms=None, onlyCA=False):
    '''

    Function to get the best superimposition between two files

    Keyword input:
    referenceFile --
    compareFile --
    mode -- 'PRECOMPUTED', 'PDBSTRING', 'PDBSTRINGBM', 'PDBSTRINGBM_RESIDUES_CONSERVED'
    DicParameters --
    algorithm -- 'nigels-core2',
    backbone --
    getDictioCorresp --
    superpose_exclude --
    listmodel --
    allAtomsList --
    allAtomsModel --
    n_iter --
    pos --
    clusterCYS --
    minmaxrms --
    onlyCA --

    Returns:
    rmsd -- rmsd between referenceFile and compareFile
    nref --
    ncom --
    allAtoms --
    compStru --
    pda --
    dictiocorres -- (only if set)

    '''

    # in mode 'PDBSTRINGBM_RESIDUES_CONSERVED', check if one of the models is larger than the other and act accordingly
    if mode == 'PDBSTRINGBM_RESIDUES_CONSERVED':
        if algorithm != "nigels-core2":
            print('The algorithm ', algorithm, ' is not supported with the mode ', mode)
            sys.exit(1)
        struref = get_structure('ref', io.StringIO(str(referenceFile)))
        strucmp = get_structure('cmp', io.StringIO(str(compareFile)))
        # file_pepe=open('ref_file.pdb','w')
        # file_pepe.write(referenceFile)
        # del file_pepe
        # file_pepa=open('comp_file.pdb','w')
        # file_pepa.write(compareFile)
        # del file_pepa
        list_atoms_ref = Bio.PDB.Selection.unfold_entities(struref, 'A')
        list_atoms_cmp = Bio.PDB.Selection.unfold_entities(strucmp, 'A')
        different_size = False
        if len(list_atoms_ref) != len(list_atoms_cmp):
            different_size = True
            # print ' \n \n ********** We entered the condition in which the two files to compare are DIFFERENT in size \n\n ***********'
            list_atoms_ref = sorted(list_atoms_ref, key=lambda x: x.get_parent().get_full_id()[3][1:])
            list_atoms_cmp = sorted(list_atoms_cmp, key=lambda x: x.get_parent().get_full_id()[3][1:])
            set_ref = set([ide.get_full_id()[3] for ide in list_atoms_ref])
            set_cmp = set([ide.get_full_id()[3] for ide in list_atoms_cmp])
            common_set = set_ref.intersection(set_cmp)
            # print 'common_set',common_set
            common_list_id = sorted(list(common_set), key=lambda x: x[1:])
            # print 'common_list_id',common_list_id
            if set_ref.issuperset(set_cmp):
                # print 'larger set corresponds to the reference. We need to modify referenceFile to the common set of atoms'
                common_list = [ele for ele in list_atoms_ref if ele.get_full_id()[3] in common_list_id]
                # print 'len(common_list)',len(common_list)
                originalreferenceFile = copy.deepcopy(referenceFile)
                referenceFile = getPDBFromListOfAtom(common_list)[0]
            else:
                # print 'larger set corresponds to the compared. We need to modify compareFile to the common set of atoms'
                common_list = [ele for ele in list_atoms_cmp if ele.get_full_id()[3] in common_list_id]
                # print 'len(common_list)',len(common_list)
                originalcompareFile = copy.deepcopy(compareFile)
                compareFile = getPDBFromListOfAtom(common_list)[0]

    if isinstance(referenceFile, str) or mode == "PRECOMPUTED":
        listValiRef = getListCA("referencia", referenceFile, mode, DicParameters=DicParameters, backbone=backbone,
                                listmodel=listmodel, algorithm=algorithm, superpose_exclude=superpose_exclude, pos=pos,
                                clusterCYS=clusterCYS)
    elif isinstance(referenceFile, list):
        listValiRef = referenceFile

    if isinstance(compareFile, str) or mode == "PRECOMPUTED":
        listValiCmp = getListCA("compare", compareFile, mode, DicParameters=DicParameters, backbone=backbone,
                                listmodel=listmodel, algorithm=algorithm, superpose_exclude=superpose_exclude, pos=pos,
                                clusterCYS=clusterCYS)
    elif isinstance(compareFile, list):
        listValiCmp = compareFile

    # print listValiRef
    # print "======================================"
    # print listValiCmp

    #  print "-----",listValiCmp
    reference = (listValiRef[0])[1]

    nref = len(reference)
    # print "algorithm",algorithm
    pdball = ""
    allAtoms = []
    compStru = None
    dictiocorresp = {}

    # print "len(listValiRef)",len(listValiRef),"len(listValiCmp)",len(listValiCmp)

    if algorithm.startswith("pymol") or algorithm.startswith("minrms") or algorithm.startswith("superpose"):
        nameficheref = "./ref_" + str(datetime.datetime.now()) + ".pdb"
        writePDBFromListOfAtom(reference, nameficheref)

        namefichecmp = "./cmp_" + str(datetime.datetime.now()) + ".pdb"
        compare = (listValiCmp[0])[1]
        writePDBFromListOfAtom(compare, namefichecmp)

        ncom = len(compare)
        # print "nref",nref,"ncom",ncom

    try:
        if algorithm == "nigels-core2":
            if isinstance(compareFile, str):
                compStru, allAtoms = getListAllAtoms("compa", compareFile)
            elif allAtomsList != None:
                allAtoms = allAtomsList
            elif allAtomsModel != None and isinstance(allAtomsModel, str):
                compStru, allAtoms = getListAllAtoms("compa", allAtomsModel)
            else:
                print(
                    "You should give a compareFile as PDB file or if compareFile is a a list of combinations of fragments as atoms, then allAtomsList should be configured")
                return (100, 0, 0, [], None, "")

            best_rmsd = 1000000000000
            best_super = None
            best_ncom = None
            best_super_ref = None
            best_nref = None
            best_R = None
            best_t = None
            howm = 1
            ifbreak = False
            tcnt = 1

            liran = list(range(len(listValiRef)))
            numpy.random.shuffle(liran)
            for cren in liran:
                nref, reference, lenfraref = listValiRef[cren]
                # print "===",[f.get_full_id() for f in nref],lenfraref
                refi = []
                k = []
                if clusterCYS:
                    reference = nref

                for atom in reference:
                    if atom.get_name() == "CA":
                        k.append(atom.get_full_id()[3][1])
                    if onlyCA and atom.get_name() != "CA":
                        continue
                    refi.append(atom.get_coord())
                refi = numpy.array(refi)

                for tren in range(len(listValiCmp)):
                    no_cmp, compare, lenfracmp = listValiCmp[tren]
                    # print "===", [f.get_full_id() for f in no_cmp], lenfracmp

                    tocontinue = False
                    if clusterCYS:
                        compare = no_cmp
                    es = list(zip(lenfraref, lenfracmp))
                    for z in range(len(es)):
                        # print "len lenfraref",len(lenfraref),"len lenfracmp",len(lenfracmp),z
                        if es[z][0] != es[z][1]:
                            tocontinue = True
                            break
                            # else:
                            #       print "--ref--",lenfraref[z],"--cmp--",lenfracmp[z]
                    # print "==================================================="
                    if tocontinue:
                        continue

                    howm += 1
                    compi = []
                    z = []
                    for atom in compare:
                        if atom.get_name() == "CA":
                            z.append(atom.get_full_id()[3][1])
                        if onlyCA and atom.get_name() != "CA":
                            continue
                        compi.append(atom.get_coord())
                    compi = numpy.array(compi)

                    transf, rmsd_list = csb.bio.utils.fit_wellordered(refi, compi, n_iter=n_iter, full_output=True, n_stdv=2,
                                                              tol_rmsd=0.005, tol_stdv=0.0005)

                    if len(rmsd_list) > 0:
                        rmsd = rmsd_list[-1][1]
                    else:
                        rmsd = 100

                    R, t = transf
                    # print "---1---============================================================"
                    # print k
                    # print "---2---==========================================================="
                    # print z
                    # print "len(refi)",len(refi),"len(compi)",len(compi),rmsd
                    ncom = len(compi)

                    if rmsd <= best_rmsd:
                        # print "rmsd",rmsd,"best_rmsd",best_rmsd
                        # print "rmsd",rmsd,"list",rmsd_list
                        outliers = 0
                        for tup in rmsd_list:
                            outliers += len(tup[2])
                            # print "==**==**==",ncom-outliers
                        tcnt = 0
                        best_rmsd = rmsd
                        best_super = tren
                        best_super_ref = cren
                        best_ncom = ncom - outliers
                        best_R = R
                        best_t = t
                    else:
                        tcnt += 1

                    if best_rmsd <= 1.0 or tcnt > 400:
                        ifbreak = True
                        break
                if ifbreak:
                    break

            # print "Comparisons done:",howm
            no_cmp, compare, lencmp = listValiCmp[best_super]
            nref, reference, lenref = listValiRef[best_super_ref]
            dictiocorresp = {}
            rmsd = best_rmsd
            ncom = best_ncom
            nref = best_ncom
            for l in range(len(reference)):
                atomr = reference[l]
                atomc = compare[l]
                rline = get_atom_line_easy(atomr)
                cline = get_atom_line_easy(atomc)
                dictiocorresp[cline[17:27]] = rline[17:27]
            if mode == 'PDBSTRINGBM_RESIDUES_CONSERVED':
                if different_size:
                    compareFile = copy.deepcopy(originalcompareFile)
                cmpStru, allAtoms = getListAllAtoms("cmpStru", compareFile)
                allAtoms = transform_atoms(allAtoms, best_R, best_t)
                pdball = getPDBFromListOfAtom(allAtoms)
                # print 'pdball',pdball
            else:
                allAtoms = transform_atoms(allAtoms, best_R, best_t)
                pdball = getPDBFromListOfAtom(allAtoms)
        elif algorithm == "nigels-core":
            if isinstance(compareFile, str):
                compStru, allAtoms = getListAllAtoms("compa", compareFile)
            elif allAtomsList != None:
                allAtoms = allAtomsList
            elif allAtomsModel != None and isinstance(allAtomsModel, str):
                compStru, allAtoms = getListAllAtoms("compa", allAtomsModel)
            else:
                print(
                    "You should give a compareFile as PDB file or if compareFile is a alist of combinations of fragments as atoms, then allAtomsList should be configured")
                return (100, 0, 0, [], None, "")

            best_rmsd = 1000000000000
            best_R = None
            best_t = None
            best_super = None
            best_ncom = None
            refi = []
            for atom in reference:
                refi.append(atom.get_coord())
            refi = numpy.array(refi)

            for tren in range(len(listValiCmp)):
                no_cmp, compare = listValiCmp[tren]
                ncom = len(compare)
                compi = []
                for atom in compare:
                    compi.append(atom.get_coord())
                compi = numpy.array(compi)

                transf, rmsd_list = csb.bio.utils.fit_wellordered(refi, compi, n_iter=n_iter, full_output=True)

                if len(rmsd_list) > 0:
                    rmsd = rmsd_list[-1][1]
                else:
                    rmsd = 100

                R, t = transf
                # print rmsd_list
                # print "================================",rmsd
                if rmsd <= best_rmsd:
                    outliers = 0
                    for tup in rmsd_list:
                        outliers += len(tup[2])
                    best_rmsd = rmsd
                    best_R = R
                    best_t = t
                    best_super = tren
                    best_ncom = ncom - outliers

            no_cmp, compare = listValiCmp[best_super]
            ncom = best_ncom
            nref = ncom
            dictiocorresp = {}
            for l in range(len(reference)):
                atomr = reference[l]
                atomc = compare[l]
                rline = get_atom_line_easy(atomr)
                cline = get_atom_line_easy(atomc)
                dictiocorresp[cline[17:27]] = rline[17:27]

            # super_imposer = Superimposer()
            # super_imposer.rotran = (best_R,best_t)
            # super_imposer.apply(allAtoms)
            allAtoms = transform_atoms(allAtoms, best_R, best_t)
            pdball = getPDBFromListOfAtom(allAtoms)
        elif algorithm == "biopython":
            if isinstance(compareFile, str):
                compStru, allAtoms = getListAllAtoms("compa", compareFile)
            elif allAtomsList != None:
                allAtoms = allAtomsList
            elif allAtomsModel != None and isinstance(allAtomsModel, str):
                compStru, allAtoms = getListAllAtoms("compa", allAtomsModel)
            else:
                print(
                    "You should give a compareFile as PDB file or if compareFile is a alist of combinations of fragments as atoms, then allAtomsList should be configured")
                return (100, 0, 0, [], None, "")

            best_rmsd = 1000000000000
            best_super = -1
            best_ncom = -1
            for tren in range(len(listValiCmp)):
                no_cmp, compare, size = listValiRef[
                    tren]  # NOTE: I am modifying this line eliminating no_cmp it might affect other functions?
                ncom = len(compare)
                # print "nref",len(reference),"ncom",len(compare)
                super_imposer = Bio.PDB.Superimposer()
                super_imposer.set_atoms(reference, compare)
                rmsd = super_imposer.rms
                # print "rmsd",rmsd,"best_rmsd",best_rmsd
                if rmsd <= best_rmsd:
                    best_rmsd = rmsd
                    best_super = tren
                    best_ncom = ncom

            no_comp, compare, size = listValiCmp[best_super]
            ncom = len(compare)
            super_imposer = Bio.PDB.Superimposer()
            super_imposer.set_atoms(reference, compare)
            dictiocorresp = {}
            for l in range(len(reference)):
                atomr = reference[l]
                atomc = compare[l]
                rline = get_atom_line_easy(atomr)
                cline = get_atom_line_easy(atomc)
                dictiocorresp[cline[17:27]] = rline[17:27]

            rmsd = super_imposer.rms
            # print "nref",len(reference),"ncom",len(compare),"rmsd",rmsd
            ncom = best_ncom
            super_imposer.apply(allAtoms)
            pdball = getPDBFromListOfAtom(allAtoms)
        elif algorithm == "biopython-core":
            if isinstance(compareFile, str):
                compStru, allAtoms = getListAllAtoms("compa", compareFile)
            elif allAtomsList != None:
                allAtoms = allAtomsList
            elif allAtomsModel != None and isinstance(allAtomsModel, str):
                compStru, allAtoms = getListAllAtoms("compa", allAtomsModel)
            else:
                print(
                    "You should give a compareFile as PDB file or if compareFile is a alist of combinations of fragments as atoms, then allAtomsList should be configured")
                return (100, 0, 0, [], None, "")

            best_rmsd = 1000000000000
            best_super = None
            best_ncom = None
            best_super_ref = None
            best_nref = None
            for cren in range(len(listValiRef)):
                nref, reference, lenfraref = listValiCmp[cren]
                for tren in range(len(listValiCmp)):
                    no_cmp, compare, lenfracmp = listValiCmp[tren]
                    tocontinue = False
                    for z in range(len(lenfraref)):
                        if lenfraref[z] != lenfracmp[z]:
                            tocontinue = True
                            break
                            # else:
                            # if lenfraref[z] not in [32,27]:
                            #       tocontinue = True
                            #       break
                            # print "--ref--",lenfraref[z],"--cmp--",lenfracmp[z]
                    # print "==================================================="
                    if tocontinue:
                        continue

                    ncom = len(compare)
                    # print "nref",len(reference),"ncom",len(compare)
                    super_imposer = Bio.PDB.Superimposer()
                    super_imposer.set_atoms(reference, compare)
                    rmsd = super_imposer.rms
                    # print "rmsd",rmsd,"best_rmsd",best_rmsd
                    if rmsd <= best_rmsd:
                        best_rmsd = rmsd
                        best_super = tren
                        best_super_ref = cren
                        best_ncom = ncom
                        best_nref = nref

            no_cmp, compare = listValiCmp[best_super]
            ncom = len(compare)
            nref, reference = listValiRef[best_super_ref]
            super_imposer = Bio.PDB.Superimposer()
            super_imposer.set_atoms(reference, compare)
            dictiocorresp = {}
            for l in range(len(reference)):
                atomr = reference[l]
                atomc = compare[l]
                rline = get_atom_line_easy(atomr)
                cline = get_atom_line_easy(atomc)
                dictiocorresp[cline[17:27]] = rline[17:27]

            rmsd = super_imposer.rms
            # print "nref",len(reference),"ncom",len(compare),"rmsd",rmsd
            ncom = best_ncom
            super_imposer.apply(allAtoms)
            pdball = getPDBFromListOfAtom(allAtoms)
        elif algorithm == "pymol_align":
            namefichepml = "./pml_" + str(datetime.datetime.now()) + ".pml"
            nameficheout = "./out_" + str(datetime.datetime.now()) + ".pdb"
            flu = open(namefichepml, "w")
            flu.write("load " + nameficheref + ", \"ref\"\n")
            flu.write("load " + namefichecmp + ", \"cmp\"\n")
            flu.write("align cmp and name ca, ref and name ca\n")
            flu.write("save " + nameficheout + ", \"cmp\"\n")
            flu.close()
            p = subprocess.Popen(['pymol', "-c", namefichepml], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            listout = out.split("\n")
            rmsd = 100
            for line in listout:
                listline = line.split()
                if line.startswith(" Executive: RMS ="):
                    rmsd = float(listline[3])
            os.remove(namefichepml)
            compStru, allAtoms = getListAllAtoms("compa", nameficheout)
            t = open(nameficheout, "r")
            pdball = t.read()
            t.close()
            os.remove(nameficheout)
        elif algorithm == "superpose":
            p = subprocess.Popen(['superpose', nameficheref, namefichecmp], stdout=subprocess.PIPE,
                                 stderr=subprocess.PIPE)
            out, err = p.communicate()
            # p2 = subprocess.Popen('superpose '+nameficheref+" "+namefichecmp+ ' | grep  "| . .\."',  shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            # out2, err2 = p2.communicate()
            rmsd = 100
            ncom = 100
            lines = out.decode().splitlines()
            for line in lines:
                if line.startswith("     r.m.s.d:"):
                    linel = line.split()
                    rmsd = float(linel[1])
                if line.startswith("      Nalign:"):
                    linel = line.split()
                    ncom = int(linel[1])
                    break

            stru_refe = get_structure("refe", nameficheref)
            stru_cmp = get_structure("cmp", namefichecmp)
            os.remove(nameficheref)
            os.remove(namefichecmp)

            dictpairing = {}
            dictpairing["reference"] = stru_refe
            dictpairing["compare"] = stru_cmp
            for line in lines:
                if len(line.split("|")) == 5 and line.split("|")[2][4] == ".":
                    chain_query = line.split("|")[1].split(":")[0][-1]
                    resid_query = int(line.split("|")[1].split(":")[1][3:].strip())
                    res_query = get_residue(stru_refe, None, chain_query, resid_query)
                    chain_target = line.split("|")[3].split(":")[0][-1]
                    resid_target = int(line.split("|")[3].split(":")[1][3:].strip())
                    res_target = get_residue(stru_cmp, None, chain_target, resid_target)
                    distance = float(line.split("|")[2][3:7])
                    dictpairing[(resid_query, chain_target, resid_target)] = (res_query, res_target, distance)
            pdball = dictpairing
        elif algorithm == "minrms":
            if minmaxrms != None:
                p = subprocess.Popen(
                    ['minrms', '-HS', '-minN', str(minmaxrms[0]), '-maxN', str(minmaxrms[1]), nameficheref,
                     namefichecmp], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            else:
                p = subprocess.Popen(['minrms', '-HS', nameficheref, namefichecmp], stdout=subprocess.PIPE,
                                     stderr=subprocess.PIPE)

            out, err = p.communicate()
            rmsd = 100
            ncom = 100
            f = open("./align_chimera.plot", "r")
            lines = f.readlines()[1:]
            f.close()
            rmsd_prev = 10000000000000
            rmsd_current = 10000000
            skiptest = False
            for line in lines:
                linel = line.split()
                os.remove("./align" + str(linel[0].strip()) + ".msf")
                if not skiptest:
                    rmsd_current = float(linel[1])
                    if rmsd_prev - rmsd_current > 0.1:
                        rmsd_prev = rmsd_current
                        ncom = int(linel[0].strip())
                    else:
                        rmsd = rmsd_prev
                        skiptest = True

            os.remove("./align_chimera.info")
            os.remove("./align_chimera.plot")
            os.remove(nameficheref)
            os.remove(namefichecmp)
            pdball = ""
        elif algorithm == "pymol_cealign":
            namefichepml = "./pml_" + str(datetime.datetime.now()) + ".pml"
            nameficheout = "./out_" + str(datetime.datetime.now()) + ".pdb"
            flu = open(namefichepml, "w")
            flu.write("load " + nameficheref + ", \"ref\"\n")
            flu.write("load " + namefichecmp + ", \"cmp\"\n")
            flu.write("cealign ref and name ca, cmp and name ca\n")
            flu.write("save " + nameficheout + ", \"cmp\"\n")
            flu.close()
            p = subprocess.Popen(['pymol', "-c", namefichepml], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            listout = out.split("\n")
            rmsd = 100
            ncom = 100
            for line in listout:
                listline = line.split()
                if line.startswith("RMSD"):
                    rmsd = float(listline[1])
                    ncom = int(listline[3])
            os.remove(namefichepml)
            compStru, allAtoms = getListAllAtoms("compa", nameficheout)
            t = open(nameficheout, "r")
            pdball = t.read()
            t.close()
            os.remove(nameficheout)
        elif algorithm == "pymol_super":
            namefichepml = "./pml_" + str(datetime.datetime.now()) + ".pml"
            nameficheout = "./out_" + str(datetime.datetime.now()) + ".pdb"
            flu = open(namefichepml, "w")
            flu.write("load " + nameficheref + ", \"ref\"\n")
            flu.write("load " + namefichecmp + ", \"cmp\"\n")
            flu.write("super cmp and name ca, ref and name ca\n")
            flu.write("save " + nameficheout + ", \"cmp\"\n")
            flu.close()
            p = subprocess.Popen(['pymol', "-c", namefichepml], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            out, err = p.communicate()
            # print "---------------------------OUT--------------------"
            # print out
            # print "---------------------------ERR--------------------"
            # print err
            # print "--------------------------------------------------"
            listout = out.split("\n")
            quit()
            rmsd = 100
            ncom = 100
            for line in listout:
                listline = line.split()
                if line.startswith(" Executive: RMS ="):
                    rmsd = float(listline[3])
                if line.startswith(" ExecutiveAlign"):
                    ncom = int(listline[1])
            os.remove(namefichepml)
            compStru, allAtoms = getListAllAtoms("compa", nameficheout)
            t = open(nameficheout, "r")
            pdball = t.read()
            t.close()
            os.remove(nameficheout)
        if algorithm.startswith("pymol"):
            os.remove(nameficheref)
            os.remove(namefichecmp)
    except:
        rmsd = 100
        print(sys.exc_info())
        traceback.print_exc(file=sys.stdout)
        # print "Incorrect File: ",compareFile
        ncom = 0
    if not getDictioCorresp:
        return (rmsd, nref, ncom, allAtoms, compStru, pdball)
    else:
        return (rmsd, nref, ncom, allAtoms, compStru, pdball, dictiocorresp)

#TODO: This function is needed by getSuperimp but must be deleted with it
def getFragmentsListWithoutCVS(modelpdb, notafile=False):
    tempDirectory = "./temp"
    if not (os.path.exists(tempDirectory)):
        os.makedirs(tempDirectory)

    allist = []

    alls = []
    globalstru = None
    if notafile:
        alls = modelpdb.splitlines()
        globalstru = get_structure("global", io.StringIO(str(modelpdb)))
    else:
        lop = open(modelpdb, "r")
        alls = lop.readlines()
        lop.close()
        globalstru = get_structure("global", modelpdb)

    countt = 0
    scrivi = True
    lai = None
    previous = -1
    previous_chain = ""
    lineaDascri = None
    for linea in alls:
        if linea.startswith("ATOM") or linea.startswith("HETATM"):
            if scrivi:
                if lai != None:
                    lai.close()
                    structure = get_structure("referencia", tempDirectory + "/" + str(countt - 1) + ".pdb")
                    reference = []
                    for model in structure.get_list():
                        for chain in model.get_list():
                            for residue in chain.get_list():
                            # atom = residue["CA"]
                                reference.append(residue)
                    allist.append(reference)
                    os.remove(tempDirectory + "/" + str(countt - 1) + ".pdb")
                lai = open(tempDirectory + "/" + str(countt) + ".pdb", "w")
                if lineaDascri != None:
                    lai.write(lineaDascri)
                    lineaDascri = None
                countt += 1
                scrivi = False
            residuo = int((linea[22:26]).strip())
            chain_res = str(linea[21])

            if previous > 0:
                resaN = get_residue(globalstru, globalstru.get_list()[0].get_id(), chain_res, (' ', residuo, ' '))
                prevResC = get_residue(globalstru, globalstru.get_list()[0].get_id(), previous_chain,
                                      (' ', previous, ' '))
                if checkContinuity(prevResC, resaN) or (residuo == previous):
                    lai.write(linea)
                    previous = residuo
                    previous_chain = chain_res
                else:
                    scrivi = True
                    lineaDascri = linea
                    previous = residuo
                    previous_chain = chain_res
            else:
                lai.write(linea)
                previous = residuo
                previous_chain = chain_res
    if lai != None:
        lai.close()
        structure = get_structure("referencia", tempDirectory + "/" + str(countt - 1) + ".pdb")
        reference = []
        for model in structure.get_list():
            for chain in model.get_list():
                for residue in chain.get_list():
                    # atom = residue["CA"]
                    reference.append(residue)
        allist.append(reference)
        os.remove(tempDirectory + "/" + str(countt - 1) + ".pdb")

    return allist

#TODO: This function is needed by getSuperimp but must be deleted with it
def createCustomFragment(structure, idname, modelV, chainV, residueList, typeOfFrag, lengthFragment):
    sst = {}
    sst["sstype"] = "nothing"
    sst["ssDescription"] = "nothing"
    residueList = sorted(residueList, key=lambda x: x[1])
    for model in structure.get_list():
        if model.get_id() != modelV:
            continue
        for chain in model.get_list():
            coordCA = []
            coordO = []
            numerOfResidues = 0
            if chain.get_id() != chainV:
                continue
            for residue in chain.get_list():
                if residue.get_id() not in residueList:
                    continue

                ca = residue["CA"]
                numerOfResidues += 1
                # per ogni atomo di Carbonio alfa, viene salvato in coordCA quattro valori:
                # [0] pos_X
                # [1] pos_Y
                # [2] pos_Z
                # [3] id del residuo di cui fa parte
                # [4] residue name code in 3 letters
                coordCA.append(
                    [float(ca.get_coord()[0]), float(ca.get_coord()[1]), float(ca.get_coord()[2]), residue.get_id(),
                     residue.get_resname()])
                ca = residue["O"]
                coordO.append(
                    [float(ca.get_coord()[0]), float(ca.get_coord()[1]), float(ca.get_coord()[2]), residue.get_id(),
                     residue.get_resname()])

            numberOfSegments = numerOfResidues - lengthFragment + 1
            # print "number of segments is",numberOfSegments
            if numberOfSegments <= 0:
                continue
            """vectorsCA = ADT.get_matrix(numberOfSegments,3)
            vectorsO = ADT.get_matrix(numberOfSegments,3)

            for i in range(lengthFragment):
                    vectorsCA[0] = [vectorsCA[0][0]+coordCA[i][0]/float(lengthFragment), vectorsCA[0][1]+coordCA[i][1]/float(lengthFragment), vectorsCA[0][2]+coordCA[i][2]/float(lengthFragment)]
                    vectorsO[0] = [vectorsO[0][0]+coordO[i][0]/float(lengthFragment), vectorsO[0][1]+coordO[i][1]/float(lengthFragment), vectorsO[0][2]+coordO[i][2]/float(lengthFragment)]

            for i in range(1,len(vectorsCA)):
                    vectorsCA[i] = [vectorsCA[i-1][0]+(coordCA[i+lengthFragment-1][0]-coordCA[i-1][0])/float(lengthFragment), vectorsCA[i-1][1]+(coordCA[i+lengthFragment-1][1]-coordCA[i-1][1])/float(lengthFragment), vectorsCA[i-1][2]+(coordCA[i+lengthFragment-1][2]-coordCA[i-1][2])/float(lengthFragment)]
                    vectorsO[i] = [vectorsO[i-1][0]+(coordO[i+lengthFragment-1][0]-coordO[i-1][0])/float(lengthFragment), vectorsO[i-1][1]+(coordO[i+lengthFragment-1][1]-coordO[i-1][1])/float(lengthFragment), vectorsO[i-1][2]+(coordO[i+lengthFragment-1][2]-coordO[i-1][2])/float(lengthFragment)]

            vectorsH = [0.0 for _ in range(numberOfSegments)]"""

            for i in range(numberOfSegments):
                # XH = vectorsCA[i][0]-vectorsO[i][0]
                # YH = vectorsCA[i][1]-vectorsO[i][1]
                # ZH = vectorsCA[i][2]-vectorsO[i][2]

                prevRes = (" ", None, " ")
                ncontigRes = 0
                resids = []
                prev_model = None
                prev_chain = None

                for yui in range(i, i + lengthFragment):  # quindi arrivo a i+lengthFragment-1
                    resan = (coordCA[yui])[3]
                    resids.append(list(resan) + [(coordCA[yui])[4]])
                    resa = resan
                    if prevRes == (" ", None, " "):
                        ncontigRes += 1
                    else:
                        resaN = get_residue(structure, model.get_id(), chain.get_id(), resan)
                        prevResC = get_residue(structure, prev_model, prev_chain, prevRes)
                        if checkContinuity(prevResC, resaN):
                            ncontigRes += 1
                    prevRes = resa
                    prev_model = model.get_id()
                    prev_chain = chain.get_id()

                if ncontigRes != lengthFragment:
                    print("The fragment to create does not contain contigous residues!")
                    return sst

                resIDs = []
                amn = []
                for q in range(lengthFragment):
                    res = coordCA[i + q][3]
                    resIDs.append(res)
                    amn.append(coordCA[i + q][4])  # take the aa name directly from the array of coordinates

                sst = recognizeFragment(idname, model, chain, coordCA[i:i + lengthFragment],
                                        coordO[i:i + lengthFragment], resIDs, amn, None, None, sstype=typeOfFrag,
                                        distS=[])

                return sst
    return sst

#TODO: This function is needed by getSuperimp but must be deleted with it
def recognizeFragment(pdb, model, chain, CAatomCoord, OatomCoord, resIDs, amn, distCV, distDF, sstype=None, distS=[]):
    if len(distS) > 0:
        if (distS[0])[1] in ["d", "bsr", "ch", "dd"]:
            distS = distS[1:]
            distCV = distCV[1:]
            distDF = distDF[1:]
            CAatomCoord = CAatomCoord[1:]
            OatomCoord = OatomCoord[1:]
            resIDs = resIDs[1:]
            amn = amn[1:]

    if len(distS) > 0:
        if (distS[-1])[1] in ["d", "bsr", "ch", "dd"]:
            distS = distS[:-1]
            distCV = distCV[:-1]
            distDF = distDF[:-1]
            CAatomCoord = CAatomCoord[:-1]
            OatomCoord = OatomCoord[:-1]
            resIDs = resIDs[:-1]
            amn = amn[:-1]

    dizio = {}
    dizio["fragLength"] = len(resIDs)
    dizio["pdbid"] = pdb
    if not isinstance(model, int):
        dizio["model"] = model.get_id()
    else:
        dizio["model"] = model

    if not isinstance(chain, str):
        dizio["chain"] = chain.get_id()
    else:
        dizio["chain"] = chain

    dizio["resIdList"] = resIDs
    dizio["amnList"] = amn
    dizio["CAatomCoord"] = CAatomCoord
    dizio["OatomCoord"] = OatomCoord

    if len(distS) == 0 and distCV != None and len(distCV) == 0:
        dizio["sstype"] = "nothing"
        dizio["ssDescription"] = "nothing"
        return dizio

    if distCV == None or distDF == None:
        # print "distribuzioni CVs non calcolate per il frammento"
        distCV = []
        distDF = []
        nS = dizio["fragLength"] - SUFRAGLENGTH + 1
        for plo in range(nS):
            # print "plo is",plo
            xca = 0.0
            yca = 0.0
            zca = 0.0
            xo = 0.0
            yo = 0.0
            zo = 0.0
            for qlo in range(SUFRAGLENGTH):
                # print "\tqlo is",qlo
                xca += CAatomCoord[plo + qlo][0]
                yca += CAatomCoord[plo + qlo][1]
                zca += CAatomCoord[plo + qlo][2]
                xo += OatomCoord[plo + qlo][0]
                yo += OatomCoord[plo + qlo][1]
                zo += OatomCoord[plo + qlo][2]
            xca /= SUFRAGLENGTH
            yca /= SUFRAGLENGTH
            zca /= SUFRAGLENGTH
            xo /= SUFRAGLENGTH
            yo /= SUFRAGLENGTH
            zo /= SUFRAGLENGTH

            XH = xca - xo
            YH = yca - yo
            ZH = zca - zo
            cv = numpy.sqrt(XH * XH + YH * YH + ZH * ZH)

            distCV.append((plo, cv))

            if plo != 0:
                distDF.append(((plo - 1, plo), numpy.abs((distCV[-1])[1] - (distCV[-2])[1])))

    dizio["distCV"] = distCV
    dizio["distDF"] = distDF

    seq = []
    for t in dizio["amnList"]:
        seq.append(AADICMAP[t])
    dizio["sequence"] = "".join(seq)

    # compute centroid CA for the fragment
    Xca = 0.0
    Yca = 0.0
    Zca = 0.0
    for i in range(len(CAatomCoord)):
        Xca += CAatomCoord[i][0]
        Yca += CAatomCoord[i][1]
        Zca += CAatomCoord[i][2]

    if len(CAatomCoord) > 0:
        Xca /= len(CAatomCoord)
        Yca /= len(CAatomCoord)
        Zca /= len(CAatomCoord)

    dizio["centroidCA"] = (Xca, Yca, Zca)

    # compute centroid O for the fragment
    Xo = 0.0
    Yo = 0.0
    Zo = 0.0
    for i in range(len(OatomCoord)):
        Xo += OatomCoord[i][0]
        Yo += OatomCoord[i][1]
        Zo += OatomCoord[i][2]

    if len(OatomCoord) > 0:
        Xo /= len(OatomCoord)
        Yo /= len(OatomCoord)
        Zo /= len(OatomCoord)

    dizio["centroidO"] = (Xo, Yo, Zo)

    # compute CV for the fragment
    XH = Xca - Xo
    YH = Yca - Yo
    ZH = Zca - Zo

    dizio["vecLength"] = numpy.sqrt(XH * XH + YH * YH + ZH * ZH)

    if sstype != None and len(distS) == len(distCV):
        if dizio["fragLength"] > 2 and sstype == "ch":
            if (distS[0])[1] == "ch":
                distS[0] = ((distS[0])[0], "dd")

            # if (distS[1])[1] == "ch":
            #        distS[1] = ((distS[1])[0],"dd")

            if (distS[-1])[1] == "ch":
                distS[-1] = ((distS[-1])[0], "dd")

            # if (distS[-2])[1] == "ch":
            #        distS[-2] = ((distS[-2])[0],"dd")

            sstype = "ah"
            for ele in distS:
                ind = ele[0]
                tip = ele[1]
                if tip == "ch":
                    sstype = "ch"
                    break

        dizio["sstype"] = sstype
        dizio["distSS"] = distS
        # print "le due lunghezze sono uguali"
        # print distS
    elif len(distCV) > 0:
        # print "ricalcolo le distSS",len(distS),len(distCV)
        meanTAH = 2.2
        meanTBS = 1.4
        st1AH = 0.10
        st2AH = 0.15
        st1BS = 0.10
        st2BS = 0.15
        outsideAlpha = 0.35
        outsideBeta = 0.30

        distSS = []

        isAlphaHelix = False
        isBetaSheet = False
        isCurvedHelix = False
        isCurvedBeta = False
        isNothing = False

        for ind in range(len(distCV)):
            (indice, cvv) = distCV[ind]

            if ind < len(distS):
                distSS.append((indice, (distS[ind])[1]))
                continue

            if numpy.abs(meanTAH - cvv) <= st1AH:
                distSS.append((indice, "ah"))
                isAlphaHelix = True
                if isBetaSheet or isCurvedBeta:
                    isNothing = True
            elif numpy.abs(meanTBS - cvv) <= st1BS:
                distSS.append((indice, "bs"))
                isBetaSheet = True
                if isAlphaHelix or isCurvedHelix:
                    isNothing = True
            elif numpy.abs(meanTAH - cvv) <= st2AH:
                distSS.append((indice, "ahd"))
                isAlphaHelix = True
                if isBetaSheet or isCurvedBeta:
                    isNothing = True
            elif numpy.abs(meanTBS - cvv) <= st2BS:
                distSS.append((indice, "bsd"))
                isBetaSheet = True
                if isAlphaHelix or isCurvedHelix:
                    isNothing = True
            else:
                distSS.append((indice, "d"))
                if isAlphaHelix:
                    isCurvedHelix = True
                    if isBetaSheet or isCurvedBeta:
                        isNothing = True
                elif isBetaSheet:
                    isCurvedBeta = True
                    if isAlphaHelix or isCurvedHelix:
                        isNothing = True

        dizio["distSS"] = distSS

        if isNothing:
            dizio["sstype"] = "nothing"
        elif isCurvedHelix:
            dizio["sstype"] = "ch"
        elif isAlphaHelix:
            dizio["sstype"] = "ah"
        elif isCurvedBeta:
            dizio["sstype"] = "cbs"
        elif isBetaSheet:
            dizio["sstype"] = "bs"
        else:
            dizio["sstype"] = "nothing"

        if sstype != None:
            dizio["sstype"] = sstype

    if dizio["sstype"] == "nothing":
        dizio["ssDescription"] = "nothing"
    elif dizio["sstype"] == "ch":
        dizio["ssDescription"] = "curved helix"
    elif dizio["sstype"] == "ah":
        dizio["ssDescription"] = "alpha helix"
    elif dizio["sstype"] == "cbs":
        dizio["ssDescription"] = "curved beta strand"
    elif dizio["sstype"] == "bs":
        dizio["ssDescription"] = "beta strand"

    return dizio

#TODO: This function is needed by getSuperimp but must be deleted with it
def getFragmentListFromListUsingSomeAtoms(struc, listaFra, caUsed, caNoUsed):
    unaOpera = True
    # print "--"
    # printSecondaryStructureElements(listaFra)
    # print "--"
    while unaOpera:
        almenoUno = True
        unaOpera = False
        while almenoUno:
            almenoUno = False
            for doing in range(len(listaFra)):
                fra1 = listaFra[doing]
                a1 = fra1["resIdList"][0]
                a1N = get_residue(struc, fra1["model"], fra1["chain"], a1)
                # a1C = getAtom(struc,fra1["model"],fra1["chain"],a1,"C")
                a2 = fra1["resIdList"][-1]
                # a2N = getAtom(struc,fra1["model"],fra1["chain"],a2,"N")
                a2C = get_residue(struc, fra1["model"], fra1["chain"], a2)
                modela2 = fra1["model"]
                chaina2 = fra1["chain"]
                pdbida2 = fra1["pdbid"]
                sst = None
                trovato = False
                for boing in range(doing + 1, len(listaFra)):
                    fra2 = listaFra[boing]
                    a3 = fra2["resIdList"][0]
                    a3N = get_residue(struc, fra2["model"], fra2["chain"], a3)
                    # a3C = getAtom(struc,fra2["model"],fra2["chain"],a3,"C")
                    a4 = fra2["resIdList"][-1]
                    # a4N = getAtom(struc,fra2["model"],fra2["chain"],a4,"N")
                    a4C = get_residue(struc, fra2["model"], fra2["chain"], a4)
                    modela3 = fra2["model"]
                    chaina3 = fra2["chain"]
                    pdbida3 = fra2["pdbid"]

                    condizio1 = modela2 == modela3 and chaina2 == chaina3 and checkContinuity(a2C, a3N)
                    condizio2 = modela2 == modela3 and chaina2 == chaina3 and checkContinuity(a4C, a1N)
                    """
                    print "paragono"
                    print fra1["resIdList"]
                    print "con"
                    print fra2["resIdList"]
                    print "condizio1",condizio1
                    print "condizio2",condizio2
                    print modela2,modela3,chaina2,chaina3,a2[1],a3[1]+1,a2[2],a3[2]
                    """
                    if condizio1:
                        sst = createCustomFragment(struc, pdbida2, modela2, chaina2,
                                                   fra1["resIdList"] + fra2["resIdList"], fra1["sstype"],
                                                   fra1["fragLength"] + fra2["fragLength"])
                        trovato = True
                    elif condizio2:
                        sst = createCustomFragment(struc, pdbida2, modela2, chaina2,
                                                   fra2["resIdList"] + fra1["resIdList"], fra2["sstype"],
                                                   fra1["fragLength"] + fra2["fragLength"])
                        trovato = True

                    if trovato:
                        break
                # print "Trovato is in frags",trovato
                if trovato:
                    almenoUno = True
                    unaOpera = True
                    nuova = []
                    # print "doing",doing,"boing",boing
                    # print "lenlistafra",len(listaFra)
                    for q in range(len(listaFra)):
                        if q not in [doing, boing]:
                            nuova.append(listaFra[q])
                            # else:
                            #       print "salto",q
                    nuova.append(sst)
                    listaFra = copy.deepcopy(nuova)
                    break

        almenoUno = True
        while almenoUno:
            almenoUno = False
            done = []
            for pos in range(len(caNoUsed)):
                elem = caNoUsed[pos]
                a1C = elem.get_parent()
                a1N = elem.get_parent()
                a1 = elem.get_parent().get_id()
                chaina1 = elem.get_parent().get_parent().get_id()
                modela1 = elem.get_parent().get_parent().get_parent().get_id()
                trovato = False
                sst = None
                for doing in range(len(listaFra)):
                    fra = listaFra[doing]
                    a2 = fra["resIdList"][0]
                    a2N = get_residue(struc, fra["model"], fra["chain"], a2)
                    # a2C = getAtom(struc,fra1["model"],fra1["chain"],a2,"C")
                    a3 = fra["resIdList"][-1]
                    # a3N = getAtom(struc,fra2["model"],fra2["chain"],a3,"N")
                    a3C = get_residue(struc, fra["model"], fra["chain"], a3)
                    modela2 = fra["model"]
                    chaina2 = fra["chain"]
                    pdbida2 = fra["pdbid"]

                    condizio1 = modela1 == modela2 and chaina1 == chaina2 and checkContinuity(a1C, a2N)
                    condizio2 = modela1 == modela2 and chaina1 == chaina2 and checkContinuity(a3C, a1N)
                    if condizio1:
                        sst = createCustomFragment(struc, pdbida2, modela2, chaina2, fra["resIdList"] + [a1],
                                                   fra["sstype"], fra["fragLength"] + 1)
                        trovato = True
                    elif condizio2:
                        sst = createCustomFragment(struc, pdbida2, modela2, chaina2, fra["resIdList"] + [a1],
                                                   fra["sstype"], fra["fragLength"] + 1)
                        trovato = True

                    if trovato:
                        done.append(pos)
                        break
                # print "Trovato is in atoms",trovato
                if trovato:
                    almenoUno = True
                    unaOpera = True
                    nuova = []
                    # print "lenlistafra",len(listaFra)
                    for q in range(len(listaFra)):
                        if q not in [doing]:
                            nuova.append(listaFra[q])
                            # else:
                            #       print "salto",q
                    nuova.append(sst)
                    listaFra = copy.deepcopy(nuova)

            news_canouse = []
            for ind in range(len(caNoUsed)):
                if ind in done:
                    caUsed.append(caNoUsed[ind])
                else:
                    news_canouse.append(caNoUsed[ind])
            caNoUsed = news_canouse

    listFragso, dicDescriptor = orderFragmentAccordingtopologicalOrder(listaFra, struc, True, True)

    if len(caNoUsed) > 0:
        print(caNoUsed)
        print("Some residues could not be added automatically to the discovered fragments.")

    # printSecondaryStructureElements(listaFra)
    # for cds in caNoUsed:
    #   print cds.get_full_id()
    # print len(caNoUsed)
    return (struc, listFragso, dicDescriptor, caUsed, caNoUsed)

#TODO: This function is needed by getSuperimp but must be deleted with it
def getFragmentListFromPDBUsingAllAtoms(inputPdbFile, drawDistri):
    tupleResult = getFragmentListFromPDB(inputPdbFile, True, drawDistri)
    struc = tupleResult[0]
    listaFra = tupleResult[1]
    caUsed = tupleResult[3]
    caNoUsed = tupleResult[4]
    # print "EEEEEEEEE",inputPdbFile
    # print "1",printSecondaryStructureElements(listaFra)
    # print "2",caNoUsed
    # print "3",len(caNoUsed)
    return getFragmentListFromListUsingSomeAtoms(struc, listaFra, caUsed, caNoUsed)

#TODO: This function is needed by getSuperimp but must be deleted with it
def getFragmentListFromPDB(pdbf, isModel, drawDistri):
    structure = None
    idname = None

    if isinstance(pdbf, str) and os.path.exists(pdbf):
        idname = os.path.basename(pdbf)[:-4]
        structure = get_structure(idname, pdbf)
    else:
        idname = "stru1"
        try:
            idname = os.path.basename(pdbf.name)[:-4]
        except:
            pass
        structure = get_structure(idname, pdbf)

    return getFragmentListFromStructure(structure, isModel, drawDistri, idname)

#TODO: This function is needed by getSuperimp but must be deleted with it
def get_matrix(n, m):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com

    :param n:
    :param m:
    :return:
    """

    return [[0.0 for _ in range(m)] for _ in range(n)]

#TODO: This function is needed by getSuperimp but must be deleted with it
def drawDistribution(fragment, windowLength):
    if not os.path.exists("./temp/drawDistri"):
        os.makedirs("./temp/drawDistri")

    ad = ""
    ae = ""
    if str(((fragment["resIdList"])[0])[2]) != " ":
        ad = str(((fragment["resIdList"])[0])[2])
    if str(((fragment["resIdList"])[-1])[2]) != " ":
        ae = str(((fragment["resIdList"])[-1])[2])

    filename = "./temp/drawDistri/" + str(fragment["pdbid"]) + "_" + str(fragment["model"]) + "_" + str(
        fragment["chain"]) + "_" + str(((fragment["resIdList"])[0])[1]) + ad + str("-") + str(
        ((fragment["resIdList"])[-1])[1]) + ae + "_" + str(fragment["sstype"]) + "_" + str(windowLength)

    qe = open(filename + ".scri", "w")
    qe.write("set terminal png size 800,800\nset output \"" + filename + ".png\"\n")
    # qe.write("set terminal canvas size 800,800\nset output \""+filename+".html\"\nplot \""+filename+".data\" using 1:2 title \"points\" with points, \""+filename+".data\" using 1:2 title \"lines\" with lines, \""+filename+".data\" using 1:2 smooth bezier title \"bezier\" with lines, \""+filename+".data\" using 1:2 smooth csplines title \"csplines\" with lines\n")

    qo = open(filename + ".data", "w")

    cont = ((fragment["resIdList"])[0])[1]
    step = 0  # windowLength-1#windowLength/2
    # print fragment["distCV"]
    # print fragment["distSS"]
    for u in range(len(fragment["distCV"])):
        if u == 0:
            qo.write("#\tX\tY\n")
        (index, jh) = (fragment["distCV"])[u]
        qo.write(" \t" + str(cont + step) + "\t" + str(jh) + "\n")
        qe.write("set label \"" + ((fragment["distSS"])[u])[1] + "\" at " + str(cont + step) + "," + str(jh) + "\n")
        cont += 1
    # qo.write(str(fragment["distCV"])+"\n")
    # qo.write(str(fragment["distSS"])+"\n")
    # qo.write(str(fragment["resIdList"])+"\n")

    qo.close()
    qe.write(
        "plot \"" + filename + ".data\" using 1:2 title \"points\" with points, \"" + filename + ".data\" using 1:2 title \"lines\" with lines, 1.4 notitle with lines, 2.2 notitle with lines\n")
    qe.close()
    # memoria = str((os.popen("~/usr/bin/gnuplot "+filename+".scri")).read())
    try:
        p = subprocess.Popen(['gnuplot', filename + ".scri"], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out, err = p.communicate()
        # os.remove(filename+".scri")
        # os.remove(filename+".data")
    except:
        print("")

    return False

#TODO: This function is needed by getSuperimp but must be deleted with it
def getFragmentListFromStructure(structure, isModel, drawDistri, idname):
    lengthFragment = SUFRAGLENGTH
    atomUsed = []
    allAtomCA = []

    h = []

    for model in structure.get_list():
        for chain in model.get_list():
            coordCA = []
            coordO = []
            numerOfResidues = 0

            for residue in chain.get_list():
                # solo se e un amminoacido ed ha almeno un carbonio alfa e un ossigeno
                if (residue.get_resname().upper() in AAList) and (residue.has_id("CA")) and (residue.has_id("O")) and (
                        residue.has_id("C")) and (residue.has_id("N")):
                    ca = residue["CA"]

                    if ca.get_occupancy() < 0.5:
                        continue
                    if (residue["O"]).get_occupancy() < 0.5:
                        continue
                    if (residue["C"]).get_occupancy() < 0.5:
                        continue
                    if (residue["N"]).get_occupancy() < 0.5:
                        continue

                    numerOfResidues += 1

                    # per ogni atomo di Carbonio alfa, viene salvato in coordCA quattro valori:
                    # [0] pos_X
                    # [1] pos_Y
                    # [2] pos_Z
                    # [3] id del residuo di cui fa parte
                    # [4] residue name code in 3 letters
                    coordCA.append(
                        [float(ca.get_coord()[0]), float(ca.get_coord()[1]), float(ca.get_coord()[2]), residue.get_id(),
                         residue.get_resname()])
                    allAtomCA.append(ca)
                    ca = residue["O"]
                    coordO.append(
                        [float(ca.get_coord()[0]), float(ca.get_coord()[1]), float(ca.get_coord()[2]), residue.get_id(),
                         residue.get_resname()])

            numberOfSegments = numerOfResidues - lengthFragment + 1
            # print numerOfResidues,lengthFragment,numberOfSegments,len(coordCA),len(coordO)
            # sys.exit(0)
            if numberOfSegments <= 0:
                # print "\t\t\tNo enough residues available to create a fragment"
                continue
            vectorsCA = get_matrix(numberOfSegments, 3)
            vectorsO = get_matrix(numberOfSegments, 3)

            for i in range(lengthFragment):
                vectorsCA[0] = [vectorsCA[0][0] + coordCA[i][0] / float(lengthFragment),
                                vectorsCA[0][1] + coordCA[i][1] / float(lengthFragment),
                                vectorsCA[0][2] + coordCA[i][2] / float(lengthFragment)]
                vectorsO[0] = [vectorsO[0][0] + coordO[i][0] / float(lengthFragment),
                               vectorsO[0][1] + coordO[i][1] / float(lengthFragment),
                               vectorsO[0][2] + coordO[i][2] / float(lengthFragment)]

            # sto cambiando che invece di -coordCA[i][0] metto coordCA[i-1][0] perche voglio togliere il precedente
            for i in range(1, len(vectorsCA)):
                vectorsCA[i] = [vectorsCA[i - 1][0] + (coordCA[i + lengthFragment - 1][0] - coordCA[i - 1][0]) / float(
                    lengthFragment),
                                vectorsCA[i - 1][1] + (coordCA[i + lengthFragment - 1][1] - coordCA[i - 1][1]) / float(
                                    lengthFragment),
                                vectorsCA[i - 1][2] + (coordCA[i + lengthFragment - 1][2] - coordCA[i - 1][2]) / float(
                                    lengthFragment)]
                vectorsO[i] = [
                    vectorsO[i - 1][0] + (coordO[i + lengthFragment - 1][0] - coordO[i - 1][0]) / float(lengthFragment),
                    vectorsO[i - 1][1] + (coordO[i + lengthFragment - 1][1] - coordO[i - 1][1]) / float(lengthFragment),
                    vectorsO[i - 1][2] + (coordO[i + lengthFragment - 1][2] - coordO[i - 1][2]) / float(lengthFragment)]

            vectorsH = [0.0 for _ in range(numberOfSegments)]

            descPosition = "" + str(idname) + "\t" + str(model.get_id()) + "\t" + str(chain.get_id()) + "\t"

            distriCV = []
            distriDf = []
            for i in range(len(vectorsCA)):
                XH = vectorsCA[i][0] - vectorsO[i][0]
                YH = vectorsCA[i][1] - vectorsO[i][1]
                ZH = vectorsCA[i][2] - vectorsO[i][2]

                prevRes = (" ", None, " ")
                ncontigRes = 0
                resids = []
                prev_model = None
                prev_chain = None

                for yui in range(i, i + lengthFragment):  # quindi arrivo a i+lengthFragment-1
                    resan = (coordCA[yui])[3]
                    resa = resan
                    resids.append(list(resan) + [(coordCA[yui])[4]])
                    if prevRes == (" ", None, " "):
                        ncontigRes += 1
                    else:
                        resaN = get_residue(structure, model.get_id(), chain.get_id(), resan)
                        prevResC = get_residue(structure, prev_model, prev_chain, prevRes)
                        # print "WHY?",resaN,prevResC,checkContinuity(prevResC,resaN)
                        if checkContinuity(prevResC, resaN):
                            ncontigRes += 1
                    prevRes = resa
                    prev_model = model.get_id()
                    prev_chain = chain.get_id()

                if ncontigRes != lengthFragment:
                    vectorsH[i] = 100  # this value identify a not reliable measure for cv
                else:
                    vectorsH[i] = numpy.sqrt(XH * XH + YH * YH + ZH * ZH)

                # print ncontigRes,lengthFragment,i
                # blocco per il calcolo delle distribuzioni
                distriCV.append((i, vectorsH[i]))
                # if vectorsH[i] == -100:
                #       print "non sono continui a",i,coordCA[i][3]
                # print "sono",i,vectorsH[i]
                if i != 0:
                    if vectorsH[i] == 100 and vectorsH[i - 1] == 100:
                        distriDf.append(((i - 1, i), 200))
                    else:
                        distriDf.append(((i - 1, i), numpy.abs(vectorsH[i] - vectorsH[i - 1])))

            resIDs = []
            amn = []
            for q in range(len(coordCA)):
                res = coordCA[q][3]
                resIDs.append(res)
                amn.append(coordCA[q][4])  # take the aa name directly from the array of coordinates
            sstp = discoverFragments(descPosition, coordCA, coordO, resIDs, structure, amn, distriCV, distriDf,
                                     lengthFragment)
            # print resIDs
            # print distriCV
            # print distriDf
            # print coordCA
            # print coordO

            for sst in sstp:
                # print sst
                # print "-----------------------------------"
                if sst["sstype"] != "nothing":
                    if sst["fragLength"] >= 3 and sst["sstype"] != "something":
                        # sst["resIdList"] = sorted(sst["resIdList"])
                        # print sst["resIdList"]
                        h.append(sst)
                        # print sst["fragLength"],"lllllllll"
                        # print sst["distCV"]
                        if drawDistri:
                            drawDistribution(sst, lengthFragment)
                        for up in sst["resIdList"]:
                            # atomUsed.append(allAtomCA[startFrag+up])
                            atomUsed.append(structure[sst["model"]][sst["chain"]][up]["CA"])
                    elif drawDistri:
                        drawDistribution(sst, lengthFragment)
                        # print sst
    # print "////////-",h
    if isModel:
        h, dicDescriptor = orderFragmentAccordingtopologicalOrder(h, structure, False, True)

        atomAvailable = []
        for u in allAtomCA:
            if u not in atomUsed:
                atomAvailable.append(u)

        return (structure, h, dicDescriptor, atomUsed, atomAvailable)
    else:
        # printSecondaryStructureElements(h)
        # sys.exit(0)
        # h = removeEqualStructures(h) #attivare l'eliminazione di strutture uguali dopo la revisione della procedura
        return (structure, h, {}, [], [])

#TODO: Rename in get_phi_psi_list
def getPhiPsiList(listResidues):
    """Return the list of phi/psi dihedral angles."""
    ppl = []
    lng = len(listResidues)
    for i in range(0, lng):
        res = listResidues[i]
        try:
            n = res['N'].get_vector()
            ca = res['CA'].get_vector()
            c = res['C'].get_vector()
        except:
            # print "Some atoms are missing"
            # print " Phi/Psi cannot be calculated for this residue"
            ppl.append((None, None))
            res.xtra["PHI"] = None
            res.xtra["PSI"] = None
            continue
        # Phi
        if i > 0:
            rp = listResidues[i - 1]
            try:
                cp = rp['C'].get_vector()
                phi = Bio.PDB.vectors.calc_dihedral(cp, n, ca, c)
                phi = phi * (180 / numpy.pi)
            except:
                phi = None
        else:
            # No phi for residue 0!
            phi = None
        # Psi
        if i < (lng - 1):
            rn = listResidues[i + 1]
            try:
                nn = rn['N'].get_vector()
                psi = Bio.PDB.vectors.calc_dihedral(n, ca, c, nn)
                psi = psi * (180 / numpy.pi)
            except:
                psi = None
        else:
            # No psi for last residue!
            psi = None
        ppl.append((phi, psi))
        # Add Phi/Psi to xtra dict of residue
        res.xtra["PHI"] = phi
        res.xtra["PSI"] = psi
    return ppl

#TODO: Rename get_residue_ramachandran_structure
def getResidueRamachandranStructure(phi, psi, residue):
    bs1Area = [(-180, 180), (-62.5, 180), (-62.5, 172.5), (-57.5, 172.5), (-57.5, 167.5), (-52.5, 167.5),
               (-52.5, 157.5), (-47.5, 157.5), (-47.5, 147.5), (-42.5, 147.5), (-42.5, 137.5), (-37.5, 137.5),
               (-37.5, 122.5), (-42.5, 122.5), (-42.5, 117.5), (-47.5, 117.5), (-47.5, 112.5), (-57.5, 112.5),
               (-57.5, 107.5), (-62.5, 107.5), (-62.5, 102.5), (-67.5, 102.5), (-67.5, 97.5), (-72.5, 97.5),
               (-72.5, 62.5), (-77.5, 62.5), (-77.5, 52.5), (-87.5, 52.5), (-87.5, 47.5), (-92.5, 47.5), (-92.5, 52.5),
               (-97.5, 52.5), (-97.5, 67.5), (-102.5, 67.5), (-102.5, 77.5), (-107.5, 77.5), (-107.5, 82.5),
               (-112.5, 82.5), (-112.5, 72.5), (-117.5, 72.5), (-117.5, 62.5), (-122.5, 62.5), (-122.5, 52.5),
               (-127.5, 52.5), (-127.5, 47.5), (-137.5, 47.5), (-137.5, 52.5), (-142.5, 52.5), (-142.5, 57.5),
               (-147.5, 57.5), (-147.5, 67.5), (-152.5, 67.5), (-152.5, 77.5), (-147.5, 77.5), (-147.5, 87.5),
               (-152.5, 87.5), (-152.5, 97.5), (-157.5, 97.5), (-157.5, 112.5), (-162.5, 112.5), (-162.5, 122.5),
               (-167.5, 122.5), (-167.5, 132.5), (-172.5, 132.5), (-172.5, 142.5), (-180, 142.5), (-180, 180)]
    bs2Area = [(-180, 180), (-42.5, 180), (-42.5, 172.5), (-42.5, 172.5), (-37.5, 172.5), (-37.5, 167.5),
               (-32.5, 167.5), (-32.5, 157.5), (-27.5, 157.5), (-27.5, 147.5), (-22.5, 147.5), (-22.5, 127.5),
               (-17.5, 127.5), (-17.5, 112.5), (-22.5, 112.5), (-22.5, 107.5), (-27.5, 107.5), (-27.5, 102.5),
               (-32.5, 102.5), (-32.5, 97.5), (-47.5, 97.5), (-47.5, 92.5), (-52.5, 92.5), (-52.5, 72.5), (-57.5, 72.5),
               (-57.5, 52.5), (-172.5, 52.5), (-177.5, 52.5), (-177.5, 77.5), (-180, 77.5), (-180, 180)]
    rah2Area = [(-57.5, 52.5), (-57.5, 42.5), (-62.5, 42.5), (-62.5, 27.5), (-57.5, 27.5), (-57.5, 22.5), (-52.5, 22.5),
                (-52.5, 12.5), (-47.5, 12.5), (-47.5, 7.5), (-42.5, 7.5), (-42.5, 2.5), (-37.5, 2.5), (-37.5, -7.5),
                (-32.5, -7.5), (-32.5, -12.5), (-27.5, -12.5), (-27.5, -27.5), (-22.5, -27.5), (-22.5, -47.5),
                (-17.5, -47.5), (-17.5, -67.5), (-22.5, -67.5), (-22.5, -77.5), (-27.5, -77.5), (-27.5, -82.5),
                (-47.5, -82.5), (-47.5, -87.5), (-77.5, -87.5), (-77.5, -92.5), (-87.5, -92.5), (-87.5, -112.5),
                (-92.5, -112.5), (-92.5, -122.5), (-97.5, -122.5), (-97.5, -137.5), (-147.5, -137.5), (-147.5, -132.5),
                (-142.5, -132.5), (-142.5, -127.5), (-147.5, -127.5), (-147.5, -97.5), (-152.5, -97.5), (-152.5, -92.5),
                (-157.5, -92.5), (-157.5, -82.5), (-162.5, -82.5), (-162.5, -52.5), (-157.5, -52.5), (-157.5, -37.5),
                (-162.5, -37.5), (-162.5, -7.5), (-167.5, -7.5), (-167.5, 32.5), (-172.5, 32.5), (-172.5, 52.5),
                (-57.5, 52.5)]
    rah1Area = [(-127.5, 47.5), (-112.5, 47.5), (-112.5, 42.5), (-102.5, 42.5), (-102.5, 37.5), (-92.5, 37.5),
                (-92.5, 32.5), (-87.5, 32.5), (-87.5, 22.5), (-82.5, 22.5), (-82.5, 17.5), (-77.5, 17.5), (-77.5, 12.5),
                (-67.5, 12.5), (-67.5, 7.5), (-62.5, 7.5), (-62.5, 2.5), (-57.5, 2.5), (-57.5, -7.5), (-52.5, -7.5),
                (-52.5, -12.5), (-47.5, -12.5), (-47.5, -22.5), (-42.5, -22.5), (-42.5, -32.5), (-37.5, -32.5),
                (-37.5, -62.5), (-42.5, -62.5), (-42.5, -67.5), (-77.5, -67.5), (-77.5, -62.5), (-117.5, -62.5),
                (-117.5, -57.5), (-122.5, -57.5), (-122.5, -47.5), (-127.5, -47.5), (-127.5, -37.5), (-132.5, -37.5),
                (-132.5, -17.5), (-137.5, -17.5), (-137.5, 2.5), (-142.5, 2.5), (-142.5, 32.5), (-137.5, 32.5),
                (-137.5, 47.5), (-127.5, 47.5)]
    # in basso a sinistra
    other1PossibleArea1 = [(-177.5, -180), (-177.5, -177.5), (-172.5, -177.5), (-172.5, -172.5), (-167.5, -172.5),
                           (-167.5, -167.5), (-127.5, -167.5), (-127.5, -172.5), (-97.5, -172.5), (-97.5, -167.5),
                           (-77.5, -167.5), (-77.5, -172.5), (-72.5, -172.5), (-72.5, -177.5), (-67.5, -177.5),
                           (-67.5, -180), (-177.5, -180)]
    # in basso a sinistra
    other1PossibleArea2 = [(-97.5, -137.5), (-92.5, -137.5), (-92.5, -142.5), (-82.5, -142.5), (-82.5, -147.5),
                           (-72.5, -147.5), (-72.5, -152.5), (-67.5, -152.5), (-67.5, -157.5), (-62.5, -157.5),
                           (-62.5, -162.5), (-57.5, -162.5), (-57.5, -167.5), (-52.5, -167.5), (-52.5, -172.5),
                           (-47.5, -172.5), (-47.5, -177.5), (-42.5, -177.5), (-42.5, -180), (-180, -180),
                           (-180, -147.5), (-97.5, -137.5), (-92.5, -137.5), (-92.5, -142.5), (-82.5, -142.5),
                           (-82.5, -147.5), (-72.5, -147.5), (-72.5, -152.5), (-67.5, -152.5), (-67.5, -157.5),
                           (-62.5, -157.5), (-62.5, -162.5), (-57.5, -162.5), (-57.5, -167.5), (-52.5, -167.5),
                           (-52.5, -172.5), (-47.5, -172.5), (-47.5, -177.5), (-42.5, -177.5), (-42.5, -180),
                           (-180, -147.5), (-177.5, -147.5), (-167.5, -147.5), (-167.5, -142.5), (-157.5, -142.5),
                           (-157.5, -137.5), (-147.5, -137.5), (-97.5, -137.5)]
    # in basso al centro
    other2PossibleArea2 = [(72.5, -102.5), (72.5, -112.5), (77.5, -112.5), (77.5, -157.5), (72.5, -157.5), (72.5, -180),
                           (57.5, -180), (57.5, -167.5), (52.5, -167.5), (52.5, -162.5), (47.5, -162.5), (47.5, -157.5),
                           (42.5, -157.5), (42.5, -152.5), (37.5, -152.5), (37.5, -142.5), (32.5, -142.5),
                           (32.5, -107.5), (37.5, -107.5), (37.5, -102.5), (42.5, -102.5), (42.5, -97.5), (52.5, -97.5),
                           (52.5, -92.5), (62.5, -92.5), (62.5, -97.5), (67.5, -97.5), (67.5, -102.5), (72.5, -102.5)]
    # in alto al centro
    other3PossibleArea2 = [(77.5, 180), (77.5, 162.5), (82.5, 162.5), (82.5, 147.5), (72.5, 147.5), (72.5, 157.5),
                           (67.5, 157.5), (67.5, 167.5), (62.5, 167.5), (62.5, 180), (77.5, 180)]
    # in alto a destra
    other4PossibleArea2 = [(162.5, 180), (162.5, 147.5), (167.5, 147.5), (167.5, 132.5), (172.5, 132.5), (172.5, 117.5),
                           (177.5, 117.5), (177.5, 77.5), (180, 77.5), (180, 180), (162.5, 180)]
    # in basso a destra
    other5PossibleArea2 = [(162.5, -180), (162.5, -177.5), (167.5, -177.5), (167.5, -167.5), (172.5, -167.5),
                           (172.5, -157.5), (177.5, -157.5), (177.5, -147.5), (180, -147.5), (180, -180), (162.5, -180)]
    lah1Area = [(57.5, 67.5), (57.5, 62.5), (62.5, 62.5), (62.5, 57.5), (67.5, 57.5), (67.5, 47.5), (72.5, 47.5),
                (72.5, 32.5), (77.5, 32.5), (77.5, 2.5), (62.5, 2.5), (62.5, 7.5), (57.5, 7.5), (57.5, 12.5),
                (52.5, 12.5), (52.5, 22.5), (47.5, 22.5), (47.5, 27.5), (42.5, 27.5), (42.5, 37.5), (37.5, 37.5),
                (37.5, 62.5), (42.5, 62.5), (42.5, 67.5), (57.5, 67.5)]
    lah2Area = [(82.5, 57.5), (87.5, 57.5), (87.5, 42.5), (92.5, 42.5), (92.5, 22.5), (97.5, 22.5), (97.5, -17.5),
                (92.5, -17.5), (92.5, -22.5), (87.5, -22.5), (87.5, -27.5), (82.5, -27.5), (82.5, -37.5), (87.5, -37.5),
                (87.5, -47.5), (92.5, -47.5), (92.5, -57.5), (87.5, -57.5), (87.5, -67.5), (82.5, -67.5), (82.5, -72.5),
                (77.5, -72.5), (77.5, -77.5), (62.5, -77.5), (62.5, -72.5), (57.5, -72.5), (57.5, -67.5), (52.5, -67.5),
                (52.5, -37.5), (57.5, -37.5), (57.5, -27.5), (62.5, -27.5), (62.5, -22.5), (57.5, -22.5), (57.5, -12.5),
                (52.5, -12.5), (52.5, -7.5), (47.5, -7.5), (47.5, -2.5), (42.5, -2.5), (42.5, 2.5), (37.5, 2.5),
                (37.5, 12.5), (32.5, 12.5), (32.5, 22.5), (27.5, 22.5), (27.5, 32.5), (22.5, 32.5), (22.5, 47.5),
                (17.5, 47.5), (17.5, 67.5), (22.5, 67.5), (22.5, 77.5), (27.5, 77.5), (27.5, 82.5), (32.5, 82.5),
                (32.5, 87.5), (47.5, 87.5), (47.5, 92.5), (67.5, 92.5), (67.5, 87.5), (72.5, 87.5), (72.5, 82.5),
                (77.5, 82.5), (77.5, 77.5), (82.5, 77.5), (82.5, 57.5)]

    # glycina in alto a sinistra
    gly1Area1 = [(-180.0, 180.0), (-180.0, 147.5), (-172.5, 147.5), (-172.5, 142.5), (-162.5, 142.5), (-162.5, 137.5),
                 (-157.5, 137.5), (-157.5, 132.5), (-137.5, 132.5), (-137.5, 127.5), (-107.5, 127.5), (-107.5, 122.5),
                 (-82.5, 122.5),
                 (-82.5, 117.5), (-77.5, 117.5), (-77.5, 112.5), (-52.5, 112.5), (-52.5, 117.5), (-47.5, 117.5),
                 (-47.5, 122.5), (-42.5, 122.5), (-42.5, 142.5), (-47.5, 142.5), (-47.5, 152.5), (-52.5, 152.5),
                 (-52.5, 162.5), (-57.5, 162.5), (-57.5, 180.0),
                 (-180.0, 180.0)]
    # glycina centro sinistra
    gly2Area1 = [(-117.5, 37.5), (-117.5, 32.5), (-122.5, 32.5), (-122.5, 27.5), (-127.5, 27.5), (-127.5, 22.5),
                 (-132.5, 22.5), (-132.5, 17.5), (-137.5, 17.5), (-137.5, 7.5), (-132.5, 7.5), (-132.5, 2.5),
                 (-127.5, 2.5), (-127.5, -7.5), (-122.5, -7.5), (-122.5, -17.5), (-117.5, -17.5), (-117.5, -22.5),
                 (-112.5, -22.5), (-112.5, -27.5), (-102.5, -27.5), (-102.5, -32.5), (-97.5, -32.5), (-97.5, -37.5),
                 (-92.5, -37.5), (-92.5, -42.5), (-87.5, -42.5), (-87.5, -47.5), (-82.5, -47.5),
                 (-82.5, -52.5), (-77.5, -52.5), (-77.5, -57.5), (-72.5, -57.5), (-72.5, -62.5), (-42.5, -62.5),
                 (-42.5, -57.5), (-37.5, -57.5), (-37.5, -42.5), (-42.5, -42.5), (-42.5, -27.5), (-47.5, -27.5),
                 (-47.5, -17.5), (-52.5, -17.5), (-52.5, -7.5), (-57.5, -7.5), (-57.5, 2.5), (-62.5, 2.5), (-62.5, 7.5),
                 (-67.5, 7.5), (-67.5, 12.5), (-72.5, 12.5), (-72.5, 17.5), (-77.5, 17.5), (-77.5, 27.5), (-82.5, 27.5),
                 (-82.5, 32.5), (-87.5, 32.5), (-87.5, 37.5), (-117.5, 37.5)]
    # glycina basso sinistra
    gly3Area1 = [(-180.0, -180.0), (-180.0, -147.5), (-162.5, -147.5), (-162.5, -152.5), (-152.5, -152.5),
                 (-152.5, -147.5), (-147.5, -147.5), (-147.5, -152.5), (-132.5, -152.5), (-132.5, -147.5),
                 (-127.5, -147.5), (-127.5, -142.5), (-122.5, -142.5), (-122.5, -137.5), (-117.5, -137.5),
                 (-117.5, -132.5), (-112.5, -132.5), (-112.5, -127.5), (-107.5, -127.5), (-107.5, -122.5),
                 (-92.5, -122.5), (-92.5, -132.5), (-87.5, -132.5), (-87.5, -137.5), (-82.5, -137.5), (-82.5, -142.5),
                 (-77.5, -142.5), (-77.5, -152.5), (-72.5, -152.5), (-72.5, -162.5), (-67.5, -162.5), (-67.5, -172.5),
                 (-62.5, -172.5), (-62.5, -177.5), (-57.5, -177.5), (-57.5, -180.0), (-180.0, -180.0)]
    # glycina alto destra
    gly4Area1 = [(180.0, 180.0), (62.5, 180.0), (62.5, 172.5), (67.5, 172.5), (67.5, 162.5), (72.5, 162.5),
                 (72.5, 152.5), (77.5, 152.5), (77.5, 142.5), (82.5, 142.5), (82.5, 137.5), (87.5, 137.5),
                 (87.5, 132.5), (92.5, 132.5), (92.5, 122.5), (107.5, 122.5), (107.5, 127.5), (112.5, 127.5),
                 (112.5, 132.5), (117.5, 132.5), (117.5, 137.5), (122.5, 137.5), (122.5, 142.5), (127.5, 142.5),
                 (127.5, 147.5), (132.5, 147.5), (132.5, 152.5), (147.5, 152.5), (147.5, 147.5), (152.5, 147.5),
                 (152.5, 152.5), (162.5, 152.5), (162.5, 147.5), (180.0, 147.5), (180.0, 180.0)]
    # glycina centro destra
    gly5Area1 = [(42.5, 62.5), (42.5, 57.5), (37.5, 57.5), (37.5, 42.5), (42.5, 42.5), (42.5, 27.5), (47.5, 27.5),
                 (47.5, 17.5), (52.5, 17.5), (52.5, 7.5), (57.5, 7.5), (57.5, -2.5), (62.5, -2.5), (62.5, -7.5),
                 (67.5, -7.5), (67.5, -12.5), (72.5, -12.5), (72.5, -17.5), (77.5, -17.5), (77.5, -27.5), (82.5, -27.5),
                 (82.5, -32.5), (87.5, -32.5), (87.5, -37.5), (117.5, -37.5), (117.5, -32.5), (122.5, -32.5),
                 (122.5, -27.5), (127.5, -27.5), (127.5, -22.5), (132.5, -22.5),
                 (132.5, -17.5), (37.5, -17.5), (137.5, -7.5), (132.5, -7.5), (132.5, -2.5), (127.5, -2.5),
                 (127.5, 7.5), (122.5, 7.5), (122.5, 17.5), (117.5, 17.5), (117.5, 22.5), (112.5, 22.5), (112.5, 27.5),
                 (102.5, 27.5), (102.5, 32.5), (97.5, 32.5),
                 (97.5, 37.5), (92.5, 37.5), (92.5, 42.5), (87.5, 42.5), (87.5, 47.5), (82.5, 47.5), (82.5, 52.5),
                 (77.5, 52.5), (77.5, 57.5), (72.5, 57.5), (72.5, 62.5), (42.5, 62.5)]
    # glycina basso centro
    gly6Area1 = [(180.0, -180.0), (57.5, -180.0), (57.5, -162.5), (52.5, -162.5), (52.5, -152.5), (47.5, -152.5),
                 (47.5, -142.5), (42.5, -142.5), (42.5, -122.5), (47.5, -122.5), (47.5, -117.5), (52.5, -117.5),
                 (52.5, -112.5), (77.5, -112.5), (77.5, -117.5), (82.5, -117.5), (82.5, -122.5), (107.5, -122.5),
                 (107.5, -127.5), (137.5, -127.5), (137.5, -132.5), (157.5, -132.5), (157.5, -137.5), (162.5, -137.5),
                 (162.5, -142.5), (172.5, -142.5), (172.5, -147.5), (180.0, -147.5), (180.0, -180.0)]
    # glycina sinistra
    gly1Area2 = [(-180.0, 180.0), (-180.0, 112.5), (-172.5, 112.5), (-172.5, 117.5), (-152.5, 117.5), (-152.5, 112.5),
                 (-137.5, 112.5), (-137.5, 107.5), (-127.5, 107.5), (-127.5, 102.5), (-117.5, 102.5), (-117.5, 97.5),
                 (-107.5, 97.5), (-107.5, 92.5), (-102.5, 92.5), (-102.5, 87.5), (-97.5, 87.5), (-97.5, 72.5),
                 (-102.5, 72.5), (-102.5, 62.5), (-107.5, 62.5), (-107.5, 57.5), (-112.5, 57.5), (-112.5, 52.5),
                 (-122.5, 52.5), (-122.5, 47.5), (-137.5, 47.5),
                 (-137.5, 52.5), (-152.5, 52.5), (-152.5, 42.5), (-157.5, 42.5), (-157.5, 7.5), (-152.5, 7.5),
                 (-152.5, -2.5), (-147.5, -2.5), (-147.5, -12.5), (-142.5, -12.5), (-142.5, -22.5), (-137.5, -22.5),
                 (-137.5, -32.5), (-132.5, -32.5),
                 (-132.5, -47.5), (-127.5, -47.5), (-127.5, -52.5), (-122.5, -52.5), (-122.5, -62.5), (-117.5, -62.5),
                 (-117.5, -72.5), (-122.5, -72.5), (-122.5, -82.5), (-127.5, -82.5), (-127.5, -87.5), (-132.5, -87.5),
                 (-132.5, -97.5), (-137.5, -97.5), (-137.5, -112.5), (-142.5, -112.5), (-142.5, -117.5),
                 (-147.5, -117.5), (-147.5, -122.5), (-152.5, -122.5), (-152.5, -127.5), (-162.5, -127.5),
                 (-162.5, -122.5), (-167.5, -122.5), (-167.5, -117.5), (-172.5, -117.5), (-172.5, -112.5),
                 (-180.0, -112.5), (-180.0, -180.0), (-47.5, -180.0), (-47.5, -167.5), (-52.5, -167.5), (-52.5, -162.5),
                 (-57.5, -162.5), (-57.5, -152.5), (-62.5, -152.5), (-62.5, -142.5), (-67.5, -142.5), (-67.5, -132.5),
                 (-72.5, -132.5), (-72.5, -127.5), (-77.5, -127.5), (-77.5, -77.5), (-42.5, -77.5), (-42.5, -72.5),
                 (-32.5, -72.5), (-32.5, -67.5), (-27.5, -67.5), (-27.5, -32.5), (-32.5, -32.5), (-32.5, -17.5),
                 (-37.5, -17.5), (-37.5, -7.5), (-42.5, -7.5), (-42.5, 2.5), (-47.5, 2.5), (-47.5, 7.5), (-52.5, 7.5),
                 (-52.5, 17.5), (-57.5, 17.5), (-57.5, 22.5), (-62.5, 22.5), (-62.5, 32.5), (-67.5, 32.5),
                 (-67.5, 72.5), (-72.5, 72.5), (-72.5, 92.5), (-67.5, 92.5), (-67.5, 97.5), (-52.5, 97.5),
                 (-52.5, 102.5), (-42.5, 102.5), (-42.5, 107.5), (-37.5, 107.5), (-37.5, 112.5), (-32.5, 112.5),
                 (-32.5, 117.5), (-27.5, 117.5), (-27.5, 142.5), (-32.5, 142.5), (-32.5, 152.5), (-37.5, 152.5),
                 (-37.5, 162.5), (-42.5, 162.5), (-42.5, 180.0), (-180.0, 180.0)]
    # glycina destra
    gly2Area2 = [(180.0, 180.0), (47.5, 180.0), (47.5, 167.5), (52.5, 167.5), (52.5, 162.5), (57.5, 162.5),
                 (57.5, 147.5), (62.5, 147.5), (62.5, 142.5), (67.5, 142.5), (67.5, 132.5), (72.5, 132.5),
                 (72.5, 127.5), (77.5, 127.5), (77.5, 77.5), (42.5, 77.5), (42.5, 72.5), (32.5, 72.5), (32.5, 67.5),
                 (27.5, 67.5), (27.5, 32.5), (32.5, 32.5), (32.5, 17.5), (37.5, 17.5), (37.5, 7.5), (42.5, 7.5),
                 (42.5, -2.5), (47.5, -2.5), (47.5, -7.5), (52.5, -7.5), (52.5, -17.5), (57.5, -17.5), (57.5, -22.5),
                 (62.5, -22.5), (62.5, -32.5), (67.5, -32.5), (67.5, -77.5), (72.5, -77.5), (72.5, -92.5),
                 (67.5, -92.5), (67.5, -97.5), (52.5, -97.5), (52.5, -102.5), (42.5, -102.5), (42.5, -107.5),
                 (37.5, -107.5), (37.5, -112.5), (37.5, -112.5), (32.5, -112.5), (32.5, -117.5), (27.5, -117.5),
                 (27.5, -142.5), (32.5, -142.5), (32.5, -152.5), (37.5, -152.5), (37.5, -162.5), (42.5, -162.5),
                 (42.5, -177.5), (47.5, -177.5), (47.5, -180.0), (180.0, -180.0),
                 (180.0, -112.5), (172.5, -112.5), (172.5, -117.5), (152.5, -117.5), (152.5, -112.5), (137.5, -112.5),
                 (137.5, -107.5), (127.5, -107.5), (127.5, -102.5), (117.5, -102.5), (117.5, -97.5), (107.5, -97.5),
                 (107.5, -92.5), (102.5, -92.5), (102.5, -87.5), (97.5, -87.5), (97.5, -67.5), (102.5, -67.5),
                 (102.5, -62.5), (107.5, -62.5), (107.5, -57.5), (112.5, -57.5), (112.5, -52.5), (122.5, -52.5),
                 (122.5, -47.5), (137.5, -47.5), (137.5, -52.5), (152.5, -52.5), (152.5, -42.5), (157.5, -42.5),
                 (157.5, -7.5), (152.5, -7.5), (152.5, 2.5), (147.5, 2.5), (147.5, 12.5), (142.5, 12.5), (142.5, 22.5),
                 (137.5, 22.5), (137.5, 32.5), (132.5, 32.5), (132.5, 47.5), (127.5, 47.5), (127.5, 52.5),
                 (122.5, 52.5), (122.5, 77.5), (127.5, 77.5), (127.5, 87.5), (132.5, 87.5), (132.5, 97.5), (37.5, 97.5),
                 (137.5, 112.5), (142.5, 112.5), (142.5, 117.5), (147.5, 117.5), (147.5, 122.5), (152.5, 122.5),
                 (152.5, 127.5), (162.5, 127.5), (162.5, 122.5), (167.5, 122.5), (167.5, 117.5), (172.5, 117.5),
                 (172.5, 112.5), (180.0, 112.5), (180.0, 180.0)]
    # prolina alto sinistra
    pro1Area1 = [(-92.500, 180.000), (-92.500, 177.500), (-97.500, 177.500), (-97.500, 142.500), (-92.500, 142.500),
                 (-92.500, 122.500), (-87.500, 122.500), (-87.500, 117.500), (-82.500, 117.500), (-82.500, 112.500),
                 (-72.500, 112.500), (-72.500, 107.500), (-67.500, 107.500), (-67.500, 112.500), (-42.500, 112.500),
                 (-42.500, 117.500), (-37.500, 117.500), (-37.500, 147.500), (-42.500, 147.500), (-42.500, 157.500),
                 (-47.500, 157.500), (-47.500, 167.500),
                 (-52.500, 167.500), (-52.500, 172.500), (-57.500, 172.500), (-57.500, 177.500), (-62.500, 177.500),
                 (-62.500, 180.000), (-92.500, 180.000)]
    # prolina basso sinistra
    pro2Area1 = [(-62.500, -180.000), (-62.500, -177.500), (-72.500, -177.500), (-72.500, -172.500),
                 (-87.500, -172.500), (-87.500, -177.500), (-92.500, -177.500), (-92.500, -180.000),
                 (-62.500, -180.000)]
    # prolina centro sinistra
    pro3Area1 = [(-72.500, 82.500), (-82.500, 82.500), (-82.500, 77.500), (-87.500, 77.500), (-87.500, 72.500),
                 (-92.500, 72.500), (-92.500, 62.500), (-87.500, 62.500), (-87.500, 52.500), (-82.500, 52.500),
                 (-82.500, 47.500), (-77.500, 47.500), (-77.500, 52.500), (-72.500, 52.500), (-72.500, 57.500),
                 (-67.500, 57.500), (-67.500, 67.500), (-72.500, 67.500), (-72.500, 82.500)]
    # prolina basso centro sinistra
    pro4Area1 = [(-92.500, 27.500), (-97.500, 27.500), (-97.500, 7.500), (-102.500, 7.500), (-102.500, -2.500),
                 (-97.500, -2.500), (-97.500, -12.500), (-92.500, -12.500), (-92.500, -22.500), (-87.500, -22.500),
                 (-87.500, -27.500), (-82.500, -27.500), (-82.500, -32.500), (-77.500, -32.500), (-77.500, -42.500),
                 (-72.500, -42.500), (-72.500, -52.500), (-67.500, -52.500), (-67.500, -57.500), (-57.500, -57.500),
                 (-57.500, -62.500), (-42.500, -62.500), (-42.500, -57.500), (-37.500, -57.500), (-37.500, -27.500),
                 (-42.500, -27.500), (-42.500, -17.500), (-47.500, -17.500), (-47.500, -12.500), (-52.500, -12.500),
                 (-52.500, -2.500), (-62.500, -2.500), (-62.500, 2.500), (-67.500, 2.500), (-67.500, 7.500),
                 (-77.500, 7.500), (-77.500, 12.500), (-82.500, 12.500), (-82.500, 17.500), (-87.500, 17.500),
                 (-87.500, 22.500), (-92.500, 22.500), (-92.500, 27.500)]
    # prolina alto sinistra
    pro1Area2 = [(-112.500, 180.000), (-112.500, 177.500), (-112.500, 142.500), (-107.500, 142.500),
                 (-107.500, 122.500), (-102.500, 122.500), (-102.500, 112.500), (-97.500, 112.500), (-97.500, 82.500),
                 (-102.500, 82.500), (-102.500, 42.500), (-107.500, 42.500), (-107.500, 37.500), (-112.500, 37.500),
                 (-112.500, -12.500), (-107.500, -12.500), (-107.500, -22.500), (-102.500, -22.500),
                 (-102.500, -32.500), (-97.500, -32.500), (-97.500, -37.500), (-92.500, -37.500), (-92.500, -42.500),
                 (-87.500, -42.500), (-87.500, -52.500), (-82.500, -52.500), (-82.500, -57.500), (-77.500, -57.500),
                 (-77.500, -62.500), (-72.500, -62.500), (-72.500, -67.500), (-62.500, -67.500), (-62.500, -72.500),
                 (-32.500, -72.500), (-32.500, -67.500), (-27.500, -67.500), (-27.500, -62.500), (-22.500, -62.500),
                 (-22.500, -32.500), (-27.500, -32.500), (-27.500, -17.500), (-32.500, -17.500), (-32.500, -7.500),
                 (-37.500, -7.500), (-37.500, -2.500), (-42.500, -2.500), (-42.500, 2.500), (-47.500, 2.500),
                 (-47.500, 7.500), (-52.500, 7.500), (-52.500, 12.500), (-57.500, 12.500), (-57.500, 17.500),
                 (-62.500, 17.500), (-62.500, 22.500), (-67.500, 22.500), (-67.500, 37.500), (-62.500, 37.500),
                 (-62.500, 47.500), (-57.500, 47.500), (-57.500, 97.500), (-47.500, 97.500), (-47.500, 102.500),
                 (-37.500, 102.500), (-37.500, 107.500), (-32.500, 107.500), (-32.500, 112.500), (-27.500, 112.500),
                 (-27.500, 152.500), (-32.500, 152.500), (-32.500, 167.500), (-37.500, 167.500), (-37.500, 172.500),
                 (-42.500, 172.500), (-42.500, 177.500), (-42.500, 180.000), (-112.500, 180.000)]
    # prolina basso sinistra
    pro2Area2 = [(-42.500, -180.000), (-42.500, -177.500), (-47.500, -177.500), (-47.500, -172.500),
                 (-52.500, -172.500), (-52.500, -167.500), (-62.500, -167.500), (-62.500, -162.500),
                 (-67.500, -162.500), (-67.500, -157.500), (-92.500, -157.500), (-92.500, -162.500),
                 (-102.500, -162.500), (-102.500, -167.500), (-107.500, -167.500), (-107.500, -177.500),
                 (-112.500, -177.500), (-112.500, -180.000), (-42.500, -180.000)]

    ritorno = "out"
    if residue.upper() == "GLY":
        if phi == None or psi == None:
            return None
        if pointInsidePolygon(phi, psi, gly1Area2):
            ritorno = "1gly2"
        if pointInsidePolygon(phi, psi, gly2Area2):
            ritorno = "2gly2"
        if pointInsidePolygon(phi, psi, gly1Area1):
            ritorno = "1gly1"
        if pointInsidePolygon(phi, psi, gly2Area1):
            ritorno = "2gly1"
        if pointInsidePolygon(phi, psi, gly3Area1):
            ritorno = "3gly1"
        if pointInsidePolygon(phi, psi, gly4Area1):
            ritorno = "4gly1"
        if pointInsidePolygon(phi, psi, gly5Area1):
            ritorno = "5gly1"
        if pointInsidePolygon(phi, psi, gly6Area1):
            ritorno = "6gly1"

        return ritorno

    if residue.upper() == "PRO":
        if phi == None or psi == None:
            return None
        if pointInsidePolygon(phi, psi, pro1Area2):
            ritorno = "1pro2"
        if pointInsidePolygon(phi, psi, pro2Area2):
            ritorno = "2pro2"
        if pointInsidePolygon(phi, psi, pro1Area1):
            ritorno = "1pro1"
        if pointInsidePolygon(phi, psi, pro2Area1):
            ritorno = "2pro1"
        if pointInsidePolygon(phi, psi, pro3Area1):
            ritorno = "3pro1"
        if pointInsidePolygon(phi, psi, pro4Area1):
            ritorno = "4pro1"

        return ritorno

    if phi == None or psi == None:
        return None
    if pointInsidePolygon(phi, psi, bs2Area):
        ritorno = "bs2"
    if pointInsidePolygon(phi, psi, rah2Area):
        ritorno = "rah2"
    if pointInsidePolygon(phi, psi, lah2Area):
        ritorno = "lah2"
    if pointInsidePolygon(phi, psi, other1PossibleArea2):
        ritorno = "1oth2"
    if pointInsidePolygon(phi, psi, bs1Area):
        ritorno = "bs1"
    if pointInsidePolygon(phi, psi, rah1Area):
        ritorno = "rah1"
    if pointInsidePolygon(phi, psi, lah1Area):
        ritorno = "lah1"
    if pointInsidePolygon(phi, psi, other1PossibleArea1):
        ritorno = "1oth1"
    if pointInsidePolygon(phi, psi, other2PossibleArea2):
        ritorno = "2oth2"
    if pointInsidePolygon(phi, psi, other3PossibleArea2):
        ritorno = "3oth2"
    if pointInsidePolygon(phi, psi, other4PossibleArea2):
        ritorno = "4oth2"
    if pointInsidePolygon(phi, psi, other5PossibleArea2):
        ritorno = "5oth2"

    return ritorno

#TODO: Rename in point_inside_polygon
def pointInsidePolygon(x, y, poly):
    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xinters = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xinters:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside

#TODO: Not used but useful: rename to compute_ramachandran and possibly improve
def computeRamachandran(structure, model, chain, resIDs, nres):
    luc = []
    tuc = []
    resi = structure[model][chain][resIDs[nres]]

    menC = False
    for residue in (structure[model][chain]):
        resaN = resi
        prevResC = residue
        if checkContinuity(prevResC, resaN):
            luc.append(structure[model][chain][residue.get_id()])
            menC = True
            break

    luc.append(resi)

    piuC = False
    for residue in (structure[model][chain]):
        resaN = residue
        prevResC = resi
        if checkContinuity(prevResC, resaN):
            luc.append(structure[model][chain][residue.get_id()])
            piuC = True
            break

    result = None
    if menC and piuC:
        # print "luc",luc
        tuc = getPhiPsiList(luc)
        # print "tuc",tuc

        phi = (tuc[1])[0]
        psi = (tuc[1])[1]
        # print "phi",phi,"psi",psi
        result = getResidueRamachandranStructure(phi, psi, luc[1].get_resname())

    # print "ededede",luc
    # print "e",luc[1]
    # print "eew",result

    return result, menC, piuC

#TODO: This function is needed by getSuperimp but must be deleted with it
def discoverFragments(descPosition, CAatomCoord, OatomCoord, resIDs, structure, amn, distCV, distDF, window):
    li = descPosition.split("\t")
    pdb = li[0]
    model = int(li[1])
    chain = li[2]

    meanTAH = 2.2
    meanTBS = 1.4
    st1AH = 0.10
    st2AH = 0.15
    st1BS = 0.10
    st2BS = 0.15
    outsideAlpha = 0.35
    outsideBeta = 0.30
    distSS = []

    isAlphaHelix = False
    isBetaSheet = False
    isCurvedHelix = False
    forceValidation = False
    startFrag = -1
    endFrag = -1
    listReturn = []

    aggio = ""
    cont = 0
    step = window - 1
    jumpindex = -1
    for ind in range(len(distCV)):
        # print "Turno:",distCV[ind]
        if jumpindex != -1 and ind < jumpindex:
            # print "salto il residuo",cont
            distSS.append((-1, ""))
            cont += 1
            continue

        jumpindex = -1
        (indice, cvv) = distCV[ind]

        if startFrag == -1:
            startFrag = ind
            # print startFrag,endFrag

            # print model,chain,cont,step

        result, menC, piuC = computeRamachandran(structure, model, chain, resIDs, cont)

        # print "per capire00",startFrag,ind
        if numpy.abs(meanTAH - cvv) <= st1AH and result not in ["bs1", "bs2", "out"]:
            distSS.append((indice, "ah"))
            isAlphaHelix = True
            if isBetaSheet:
                endFrag = ind + step
                """
                print "creo frammento ",cont,"modo a"
                print "==============================0",startFrag,ind,ind-1
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]
                print amn[startFrag:endFrag]
                """
                diston = distSS[startFrag:ind]
                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype="bs", distS=diston)
                listReturn.append(st)
                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
        elif numpy.abs(meanTBS - cvv) <= st1BS and result not in ["rah1", "rah2", "lah1", "lah2", "out"]:
            distSS.append((indice, "bs"))
            # distCV[ind] = (indice,cvv,"bs")
            isBetaSheet = True
            if isAlphaHelix or isCurvedHelix:
                # print resIDs[startFrag:endFrag]
                # print isAlphaHelix,isCurvedHelix,isBetaSheet
                # sys.exit(0)
                ti = ""
                if isCurvedHelix:
                    ti = "ch"
                elif isAlphaHelix:
                    ti = "ah"

                endFrag = ind + step
                """print "creo frammento ",cont,"modo b"
                print "==============================0"
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]
                print amn[startFrag:endFrag]"""

                diston = distSS[startFrag:ind]

                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype=ti, distS=diston)
                listReturn.append(st)
                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
        elif numpy.abs(meanTAH - cvv) <= st2AH and result not in ["bs1", "bs2", "out"]:
            distSS.append((indice, "ahd"))
            # distCV[ind] = (indice,cvv,"ahd")
            isAlphaHelix = True
            if isBetaSheet:
                endFrag = ind + step
                """
                print "creo frammento ",cont,"modo c"
                print "==============================0"
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]
                print amn[startFrag:endFrag]
                """
                diston = distSS[startFrag:ind]

                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype="bs", distS=diston)
                listReturn.append(st)
                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
        elif numpy.abs(meanTBS - cvv) <= st2BS and result not in ["rah1", "rah2", "lah1", "lah2", "out"]:
            distSS.append((indice, "bsd"))
            # distCV[ind] = (indice,cvv,"bsd")
            isBetaSheet = True
            if isAlphaHelix or isCurvedHelix:
                endFrag = ind + step
                ti = ""
                if isCurvedHelix:
                    ti = "ch"
                elif isAlphaHelix:
                    ti = "ah"
                """print "creo frammento ",cont,"modo d"
                print "==============================0"
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]"""

                diston = distSS[startFrag:ind]

                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype=ti, distS=diston)
                listReturn.append(st)
                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
        else:
            # print "/////////"
            # print indice
            # print cvv
            # print "////////"
            distSS.append((indice, "d"))
            # print "distorsioneeeeeee"
            # print numpy.abs(meanTAH-cvv),st2AH,meanTAH,cvv
            if cvv == 100 or ((ind + step) - startFrag) + step >= 3:
                endFrag = ind + step
                # print "tagliato",cvv,startFrag,endFrag
                if endFrag != startFrag:
                    forceValidation = True
                    if isCurvedHelix:
                        aggio = "ch"
                    elif isAlphaHelix:
                        aggio = "ah"
                    else:
                        aggio = "bs"
                # print ind,aggio
                """
                elif (((isAlphaHelix or isCurvedHelix) and numpy.abs(meanTAH-cvv) >= outsideAlpha) or (isBetaSheet and numpy.abs(meanTBS-cvv) >= outsideBeta)):
                endFrag = ind+step
                forceValidation = True
                if isCurvedHelix:
                        aggio = "ch"
                elif isAlphaHelix:
                        aggio = "ah"
                else:
                        aggio = "bs"
            elif (((isAlphaHelix or isCurvedHelix) and numpy.abs(meanTAH-cvv) < outsideAlpha and ind == len(distCV)-1) or (isBetaSheet and numpy.abs(meanTBS-cvv) >= outsideBeta and ind == len(distCV)-1)):
                endFrag = ind+step
                forceValidation = True
                if isCurvedHelix:
                        aggio = "ch"
                elif isAlphaHelix:
                        aggio = "ah"
                else:
                        aggio = "bs"
            elif not isAlphaHelix and not isCurvedHelix and not isBetaSheet:
                endFrag = ind+step
                forceValidation = True
                aggio = "nothing"
                """

            toContinue = True
            resultL = []
            result = None
            menC = False
            piuC = False
            for wer in range(window):
                result, menC, piuC = computeRamachandran(structure, model, chain, resIDs, cont + wer)
                resultL.append(result)
                if result in [None, "out"]:
                    # print "=??????",result,resIDs,cont+wer,menC,piuC
                    toContinue = False
                    break

            # print "helix",isAlphaHelix,"curved",isCurvedHelix,"beta",isBetaSheet,"FV",forceValidation,"CONTINUE",toContinue
            if toContinue:
                # print "cococococococococococcccccccccccccccccccccccccccccccccccccccccccccccc"
                # print "result",result,menC,piuC
                # print "prima",startFrag,endFrag
                if forceValidation:
                    if result == None and not piuC and ind == len(distCV) - 1:
                        endFrag = ind + step + 1
                        ind += 1
                    else:
                        endFrag = ind + step
                    # print "dopo",startFrag,endFrag
                    distSS[-1] = (indice, "ch")

                    """
                    print "creo frammento ",cont,"modo l"
                    print "==============================0"
                    print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                    print distCV[startFrag:ind]
                    print distSS[startFrag:ind]
                    print resIDs[startFrag:endFrag]
                    """
                    diston = distSS[startFrag:ind]

                    st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag],
                                           OatomCoord[startFrag:endFrag], resIDs[startFrag:endFrag],
                                           amn[startFrag:endFrag], distCV[startFrag:ind], distDF[startFrag:ind],
                                           sstype=aggio, distS=diston)
                    listReturn.append(st)
                    # print st
                    # print "//////////"
                    # sys.exit(0)
                    aggio = ""
                    if st["fragLength"] >= 3:
                        startFrag = ind + step
                        endFrag = -1
                        jumpindex = ind + step
                        # print "=0=",startFrag
                    else:
                        startFrag = -1
                        endFrag = -1
                        jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
                if not forceValidation:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1
                forceValidation = False
                cont += 1
                # print "subframmento fuori ramachandran"
                # print "ind post",ind
                # print "=1=",startFrag,endFrag
                continue

            result = resultL[0]

            if (isAlphaHelix or isCurvedHelix) and result not in ["bs1", "bs2"]:  # and ramachandran dice che e ah
                isCurvedHelix = True
                distSS[-1] = (indice, "ch")
                """
                if isBetaSheet:
                    endFrag = ind+step

                    print "creo frammento ",cont,"modo e"
                    print "==============================0"
                    print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                    print distCV[startFrag:ind]
                    print distSS[startFrag:ind]
                    print resIDs[startFrag:endFrag]

                    diston = distSS[startFrag:ind]
                    st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag], resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind], distDF[startFrag:ind], sstype="bs", distS=diston)
                    listReturn.append(st)
                    startFrag = ind+step
                    endFrag = -1
                    jumpindex = ind+step
                    isAlphaHelix = False
                    isBetaSheet = False
                    isCurvedHelix = False
                    #print "trovata elica durante una beta sheet"
                """
            elif (isAlphaHelix or isCurvedHelix) and result in ["bs1", "bs2"]:
                ti = ""
                if isCurvedHelix:
                    ti = "ch"
                elif isAlphaHelix:
                    ti = "ah"

                # distSS[-1] = (indice,ti)
                endFrag = ind + step

                """print "creo frammento ",cont,"modo f"
                print "==============================0"
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]"""

                diston = distSS[startFrag:ind]

                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype=ti, distS=diston)
                listReturn.append(st)
                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
                # print "trovata una beta o un out dentro un elica"
            elif isBetaSheet and result in ["rah1", "rah2", "lah1", "lah2"]:
                endFrag = ind + step
                # print "prima",startFrag,endFrag
                # distSS[-1] = (indice,"ch")
                """print "creo frammento ",cont,"modo g"
                print "==============================0"
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]"""
                diston = distSS[startFrag:ind]

                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype="bs", distS=diston)
                listReturn.append(st)
                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False
                # print "trovata un elica en una beta strand",result,cont,step,window
                # print "dopo",startFrag,endFrag
            elif isBetaSheet and result not in ["rah1", "rah2", "lah1", "lah2"]:
                distSS[-1] = (indice, "bsr")
            elif result in ["bs1", "bs2"]:
                isBetaSheet = True
                distSS[-1] = (indice, "bsr")
            elif result in ["rah1", "rah2", "lah1", "lah2"]:
                isCurvedHelix = True
                distSS[-1] = (indice, "ch")

            if result in ["1gly2", "2gly2", "1gly1", "2gly1", "3gly1", "4gly1", "5gly1", "6gly1", "1pro2", "2pro2",
                          "1pro1", "2pro1", "3pro1", "4pro1"]:
                distSS[-1] = (indice, result)
                # distSS[-1] = (indice,(distSS[-1])[1]+"_"+result)

            if forceValidation:
                endFrag = ind + step
                # print "prima",startFrag,endFrag
                distSS[-1] = (indice, "ch")
                """print "creo frammento ",cont,"modo z"
                print "=====================lklklk=========0"
                print len(distCV[startFrag:ind]),len(distSS[startFrag:ind])
                print distCV[startFrag:ind]
                print distSS[startFrag:ind]
                print resIDs[startFrag:endFrag]"""
                diston = distSS[startFrag:ind]

                st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:endFrag], OatomCoord[startFrag:endFrag],
                                       resIDs[startFrag:endFrag], amn[startFrag:endFrag], distCV[startFrag:ind],
                                       distDF[startFrag:ind], sstype=aggio, distS=diston)
                listReturn.append(st)
                aggio = ""
                forceValidation = False
                isAlphaHelix = False
                isBetaSheet = False
                isCurvedHelix = False

                if st["fragLength"] >= 3:
                    startFrag = ind + step
                    endFrag = -1
                    jumpindex = ind + step
                else:
                    startFrag = -1
                    endFrag = -1
                    jumpindex = -1

                    # print "subframmento fuori ramachandran"
                    # print "ind post",ind
        cont += 1

    # print "fuoriiiiii",ind,cont
    if startFrag != -1 and len(distCV[startFrag:]) > 0 and (isCurvedHelix or isAlphaHelix or isBetaSheet):
        tip = ""
        if isCurvedHelix:
            tip = "ch"
        elif isAlphaHelix:
            tip = "ah"
        else:
            tip = "bs"

        """print "creo frammento ",cont,"modo h"
        print "==============================0"
        print len(distCV[startFrag:]),len(distSS[startFrag:])
        print distCV[startFrag:]
        print distSS[startFrag:]
        print resIDs[startFrag:]"""
        diston = distSS[startFrag:]

        st = recognizeFragment(pdb, model, chain, CAatomCoord[startFrag:], OatomCoord[startFrag:], resIDs[startFrag:],
                               amn[startFrag:], distCV[startFrag:], distDF[startFrag:], sstype=tip, distS=diston)
        listReturn.append(st)

    # printSecondaryStructureElements(listReturn)

    h = listReturn
    # refino and joining de las estructuras,
    More = True
    while More:
        More = False
        for a1 in range(len(h)):
            if len((h[a1])["resIdList"]) == 0:
                continue
            for a2 in range(len(h)):
                if a1 == a2 or len((h[a2])["resIdList"]) == 0:
                    continue
                a1S = ((h[a1])["resIdList"])[0]
                a1E = ((h[a1])["resIdList"])[(h[a1])["fragLength"] - 1]
                a2S = ((h[a2])["resIdList"])[0]
                a2E = ((h[a2])["resIdList"])[(h[a2])["fragLength"] - 1]

                a1SN = get_residue(structure, (h[a1])["model"], (h[a1])["chain"], a1S)
                a1EC = get_residue(structure, (h[a1])["model"], (h[a1])["chain"], a1E)
                a2SN = get_residue(structure, (h[a2])["model"], (h[a2])["chain"], a2S)
                a2EC = get_residue(structure, (h[a2])["model"], (h[a2])["chain"], a2E)

                condizio0 = ((h[a1])["model"] == (h[a2])["model"]) and ((h[a1])["chain"] == (h[a2])["chain"])
                condizio1 = checkContinuity(a1EC, a2SN)
                condizio2 = checkContinuity(a2EC, a1SN)

                condizio4 = (h[a1])["fragLength"] >= 3 and (h[a2])["fragLength"] >= 3

                if not (condizio0 and condizio4 and (condizio1 or condizio2)):
                    continue

                condizio3 = ((h[a1])["sstype"] == (h[a2])["sstype"] == "bs") or (
                    (h[a1])["sstype"] == "cbs" and (h[a2])["sstype"] == "bs") or (
                                (h[a1])["sstype"] == "bs" and (h[a2])["sstype"] == "cbs") or (
                                (h[a1])["sstype"] == (h[a2])["sstype"] == "cbs")

                condizio8 = ((h[a1])["sstype"] == (h[a2])["sstype"] == "ah") or (
                    (h[a1])["sstype"] == "ch" and (h[a2])["sstype"] == "ah") or (
                                (h[a1])["sstype"] == "ah" and (h[a2])["sstype"] == "ch") or (
                                (h[a1])["sstype"] == (h[a2])["sstype"] == "ch")

                treshold = 0
                if condizio3:
                    treshold = 30
                elif condizio8:
                    treshold = 35

                mod1cCA = (h[a1])["centroidCA"]
                mod1cO = (h[a1])["centroidO"]
                mod2cCA = (h[a2])["centroidCA"]
                mod2cO = (h[a2])["centroidO"]
                X1 = mod1cCA[0] - mod1cO[0]
                Y1 = mod1cCA[1] - mod1cO[1]
                Z1 = mod1cCA[2] - mod1cO[2]
                X2 = mod2cCA[0] - mod2cO[0]
                Y2 = mod2cCA[1] - mod2cO[1]
                Z2 = mod2cCA[2] - mod2cO[2]

                TetaReal = angle_between([X1, Y1, Z1], [X2, Y2, Z2], [1.0, 1.0, 1.0], signed=False)
                TetaDeg = TetaReal * 57.2957795
                condizio5 = TetaDeg <= treshold

                # print "Possibilita di unire eliche",condizio8,condizio5,h[a1]["resIdList"][0],h[a2]["resIdList"][0]
                if condizio3 and condizio5:
                    tipo = "cbs"  # this is useful for bs
                    # descript = "beta sheet"

                    sstoc = None
                    if condizio1:  # a1--a2
                        sstoc = recognizeFragment((h[a1])["pdbid"], (h[a1])["model"], (h[a1])["chain"],
                                                  (h[a1])["CAatomCoord"] + (h[a2])["CAatomCoord"],
                                                  (h[a1])["OatomCoord"] + (h[a2])["OatomCoord"],
                                                  (h[a1])["resIdList"] + (h[a2])["resIdList"],
                                                  (h[a1])["amnList"] + (h[a2])["amnList"], None, None, sstype=tipo,
                                                  distS=(h[a1])["distSS"] + (h[a2])["distSS"])
                    elif condizio2:  # a2--a1
                        sstoc = recognizeFragment((h[a2])["pdbid"], (h[a2])["model"], (h[a2])["chain"],
                                                  (h[a2])["CAatomCoord"] + (h[a1])["CAatomCoord"],
                                                  (h[a2])["OatomCoord"] + (h[a1])["OatomCoord"],
                                                  (h[a2])["resIdList"] + (h[a1])["resIdList"],
                                                  (h[a2])["amnList"] + (h[a1])["amnList"], None, None, sstype=tipo,
                                                  distS=(h[a2])["distSS"] + (h[a1])["distSS"])

                    h = SystemUtility.multi_delete(h, [a1, a2])
                    h.append(sstoc)
                    More = True
                    break
                elif condizio8 and condizio5:
                    tipo = "ch"
                    sstoc = None
                    if condizio1:  # a1--a2
                        sstoc = recognizeFragment((h[a1])["pdbid"], (h[a1])["model"], (h[a1])["chain"],
                                                  (h[a1])["CAatomCoord"] + (h[a2])["CAatomCoord"],
                                                  (h[a1])["OatomCoord"] + (h[a2])["OatomCoord"],
                                                  (h[a1])["resIdList"] + (h[a2])["resIdList"],
                                                  (h[a1])["amnList"] + (h[a2])["amnList"], None, None, sstype=tipo,
                                                  distS=(h[a1])["distSS"] + (h[a2])["distSS"])
                    elif condizio2:  # a2--a1
                        sstoc = recognizeFragment((h[a2])["pdbid"], (h[a2])["model"], (h[a2])["chain"],
                                                  (h[a2])["CAatomCoord"] + (h[a1])["CAatomCoord"],
                                                  (h[a2])["OatomCoord"] + (h[a1])["OatomCoord"],
                                                  (h[a2])["resIdList"] + (h[a1])["resIdList"],
                                                  (h[a2])["amnList"] + (h[a1])["amnList"], None, None, sstype=tipo,
                                                  distS=(h[a2])["distSS"] + (h[a1])["distSS"])

                    h = SystemUtility.multi_delete(h, [a1, a2])
                    h.append(sstoc)
                    More = True
                    break

            if More:
                break
    return h

#TODO: This function is needed by getSuperimp but must be deleted with it
def angle_between(A, B, N, signed=True):
    # ANGLE BETWEEN TWO 3D VECTORS:
    # 1- dot(norm(A),norm(B)) (ANGLES UNSIGNED, PROBLEMS FOR SMALL ANGLES WITH ROUNDINGS)
    # 2- arcos(dot(A,B)/(|A|*|B|))  (ANGLE UNSIGNED, PROBLEMS FOR SMALL ANGLES WITH ROUNDINGS)
    # 3- arctan2(|cross(A,B)|,dot(A,B)) (ANGLE UNSIGNED BUT NOT PROBLEMS OF ROUNDINGS
    #   define a vector NORM ex.: N = [0,0,1]
    #   sign = dot(NORM,cross(A,B))
    #   if sign < 0 then ANGLE measured in 3 should be negative

    CrossX = A[1] * B[2] - A[2] * B[1]
    CrossY = A[2] * B[0] - A[0] * B[2]
    CrossZ = A[0] * B[1] - A[1] * B[0]

    fCross = numpy.sqrt(CrossX * CrossX + CrossY * CrossY + CrossZ * CrossZ)
    scaP2 = (A[0] * B[0]) + (A[1] * B[1]) + (A[2] * B[2])
    Teta_2 = numpy.arctan2(fCross, scaP2)

    if signed:
        sign = (N[0] * CrossX) + (N[1] * CrossY) + (N[2] * CrossZ)
        if sign < 0:
            Teta_2 = -Teta_2

        return Teta_2
    else:
        return Teta_2

#TODO: This function is needed by getSuperimp but must be deleted with it
def computeDistancesForFragments(listfrag):
    dicDescriptor = {}
    maxDist = numpy.NINF
    minDist = numpy.inf
    h = listfrag

    for i in range(len(h)):
        if "distances" in (h[i]).keys():
            (h[i])["distances"] = []

    for i in range(len(h)):
        maxDf = numpy.NINF
        minDf = numpy.inf
        minIf = -1
        maxIf = -1
        maxCf = (-1, -1, -1)
        minCf = (-1, -1, -1)
        for j in range(len(h)):
            if j == i:
                continue
            mod1cCA = (h[i])["centroidCA"]
            mod2cCA = (h[j])["centroidCA"]
            X = mod1cCA[0] - mod2cCA[0]
            Y = mod1cCA[1] - mod2cCA[1]
            Z = mod1cCA[2] - mod2cCA[2]
            di = numpy.sqrt(X * X + Y * Y + Z * Z)
            if di > maxDist:
                maxDist = di

            if di < minDist:
                minDist = di

            if di > maxDf:  # and i < j:
                # if j < i and "distances" in h[j] and i == (((h[j])["distances"])[0])[0]:
                #       continue

                maxDf = di
                maxIf = j
                maxCf = (X, Y, Z)

            if di < minDf:  # and i < j:
                # if j < i and "distances" in h[j] and i == (((h[j])["distances"])[1])[0]:
                #       continue

                minDf = di
                minIf = j
                minCf = (X, Y, Z)

        if "distances" not in (h[i]).keys():
            (h[i])["distances"] = []
            # if i != len(h)-1:
            (h[i])["distances"].append([maxIf, maxCf, maxCf])  # (nmFrag,distance,direction)
            (h[i])["distances"].append([minIf, minCf, minCf])  # (nmFrag,distance,direction)
            # else:
            #       (h[i])["distances"] = []
        else:
            # if i != len(h)-1:
            (h[i])["distances"].append([maxIf, maxCf, maxCf])  # (nmFrag,distance,direction)
            (h[i])["distances"].append([minIf, minCf, minCf])  # (nmFrag,distance,direction)
            # else:
            #       (h[i])["distances"] = []

    dicDescriptor["rangeDistances"] = (minDist, maxDist)

    return h, dicDescriptor

#TODO: Remove from ADT
def strongly_connected_components(graph):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com

    :param graph:
    :return:
    """

    result = []
    stack = []
    low = {}

    def visit(node):
        if node in low: return

        num = len(low)
        low[node] = num
        stack_pos = len(stack)
        stack.append(node)

        for successor in graph[node]:
            visit(successor)
            low[node] = min(low[node], low[successor])

        if num == low[node]:
            component = tuple(stack[stack_pos:])
            del stack[stack_pos:]
            result.append(component)
            for item in component:
                low[item] = len(graph)

    for node in graph:
        visit(node)

    return result

#TODO: Remove from ADT
def topological_sort(graph):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com

    :param graph:
    :return:
    """
    count = {}
    for node in graph:
        count[node] = 0
    for node in graph:
        for successor in graph[node]:
            count[successor] += 1

    ready = [node for node in graph if count[node] == 0]

    result = []
    while ready:
        node = ready.pop(-1)
        result.append(node)

        for successor in graph[node]:
            count[successor] -= 1
            if count[successor] == 0:
                ready.append(successor)

    return result

#TODO: Remove from ADT
def robust_topological_sort(graph):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com

    :param graph:
    :return:
    """

    components = strongly_connected_components(graph)

    node_component = {}
    for component in components:
        for node in component:
            node_component[node] = component

    component_graph = {}
    for component in components:
        component_graph[component] = []

    for node in graph:
        node_c = node_component[node]
        for successor in graph[node]:
            successor_c = node_component[successor]
            if node_c != successor_c:
                component_graph[node_c].append(successor_c)

    return topological_sort(component_graph)

#TODO: This function is needed by getSuperimp but must be deleted with it
def orderFragmentAccordingtopologicalOrder(FragmentList, structure, removeEqualStructure, isModel):
    h = FragmentList
    h, dicDescriptor = computeDistancesForFragments(h)

    # printSecondaryStructureElements(h)
    for i in range(len(h)):
        tg = h[i]
    # print "frag ",i,"-----"
    #       for kk in tg["distances"]:
    #               print kk

    dizzz = {}

    for j in range(len(h)):
        liu = []
        for i in range(len(h)):
            tg = h[i]
            for lup in range(len(tg["distances"])):
                # if lup == 0: #for the ordering i consider only the minimum distance
                #       continue
                kk = (tg["distances"])[lup]
                if kk[0] == j:
                    liu.append(i)
        dizzz[j] = liu

    nuovoOrd = robust_topological_sort(dizzz)
    # print nuovoOrd

    """h2 = []
    h3 = deepcopy(h)
    cont = 0
    for u in nuovoOrd:
            for k in u:
                    #(h[k])["distances"] = []
                    for s in range(len(h)):
                            for q in range(len((h[s])["distances"])):
                                    elem = ((h[s])["distances"])[q]
                                    if elem[0] == k:
                                            (((h3[s])["distances"])[q])[0] = cont
                    h2.append(h3[k])
                    cont += 1"""

    h2 = []
    for u in nuovoOrd:
        for k in u:
            h2.append(h[k])

    h2, dicDescriptor = computeDistancesForFragments(h2)

    dizzz = {}

    for j in range(len(h2)):
        liu = []
        for i in range(len(h2)):
            tg = h2[i]
            for lup in range(len(tg["distances"])):
                # if lup == 0: #for the ordering i consider only the minimum distance
                #       continue
                kk = (tg["distances"])[lup]
                if kk[0] == j:
                    liu.append(i)
        dizzz[j] = liu

    nuovoOrd = robust_topological_sort(dizzz)
    # print nuovoOrd

    for i in range(len(h2)):
        tg = h2[i]
    # print "frag ",i,"-----"
    #       for kk in tg["distances"]:
    #               print kk
    return h2, dicDescriptor

def get_atom_line_easy(atom):
    item = atom
    orig_atom_num = item.get_serial_number()
    hetfield, resseq, icode = item.get_parent().get_id()
    segid = item.get_parent().get_segid()
    resname = item.get_parent().get_resname()
    chain_id = item.get_parent().get_parent().get_id()
    element = item.get_name()
    return get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id)

@SystemUtility.deprecated("Use get_pdb_from_list_of_frags")
def getPDBFormattedAsString(nModel, lista, strucActualPDB, pathBase, dizioConv={}, externalRes=[], useDizioConv=True,
                            normalize=False, bfactorNor=25.0):
    pdbid = str(strucActualPDB.get_id())
    nomeFilefine = pathBase + "/" + pdbid + "_"
    myHea = Heap()
    consumedResidue = []

    def actionImplemented(ModelT, ChainT, rex, consumedResidue, myHea):
        for model in strucActualPDB.get_list():
            if model.get_id() != ModelT:
                continue
            for chain in model.get_list():
                if chain.get_id() != ChainT:
                    continue
                for residue in chain.get_list():
                    if residue.get_id() != rex:
                        continue

                    """if residue.get_id() == (' ', 186, ' ') and dizioConv == {}:
                            print "residue",residue.get_id()
                            print "residue full id",residue.get_full_id()
                            print "model",model.get_id()
                            print "chain",chain.get_id()
                            print "is not in consumedResidue",residue.get_full_id() not in consumedResidue"""

                    if residue.get_full_id() not in consumedResidue:
                        if residue.has_id("N"):
                            pos = 0
                            myHea.push((model.get_id(), chain.get_id(),
                                        (residue.get_id()[0], residue.get_id()[1], residue.get_id()[2]), pos),
                                       (residue["N"], "N"))
                        if residue.has_id("CA"):
                            pos = 1
                            myHea.push((model.get_id(), chain.get_id(),
                                        (residue.get_id()[0], residue.get_id()[1], residue.get_id()[2]), pos),
                                       (residue["CA"], "CA"))
                        if residue.has_id("CB"):
                            pos = 2
                            myHea.push((model.get_id(), chain.get_id(),
                                        (residue.get_id()[0], residue.get_id()[1], residue.get_id()[2]), pos),
                                       (residue["CB"], "CB"))
                        if residue.has_id("C"):
                            pos = 3
                            myHea.push((model.get_id(), chain.get_id(),
                                        (residue.get_id()[0], residue.get_id()[1], residue.get_id()[2]), pos),
                                       (residue["C"], "C"))
                        if residue.has_id("O"):
                            pos = 4
                            myHea.push((model.get_id(), chain.get_id(),
                                        (residue.get_id()[0], residue.get_id()[1], residue.get_id()[2]), pos),
                                       (residue["O"], "O"))

                        consumedResidue.append(residue.get_full_id())
        return (consumedResidue, myHea)


        # in this way for each residue i will append only the CA atom that is undersored because
        # Biopython gives me as default the atom with major occupancy when i call residue["CA"],
        # moreover i don't worry about the  < 0.5 occupancy exlusion criteria because if the residue
        # it is included in the solution is because his CA atom has an occupancy > 0.5.

    for resi in externalRes:
        (consumedResidue, myHea) = actionImplemented(resi[0], resi[1], resi[2], consumedResidue, myHea)

    for y in range(len(lista)):
        fragment = lista[y]
        for i in range(len(fragment["resIdList"])):
            re = (fragment["resIdList"])[i]
            if i == 0:
                # nomeFilefine += fragment["chain"]+str(re)
                if len(dizioConv.keys()) > 0:
                    nomeFilefine += str((dizioConv[(fragment["chain"], re, "CA")])[1]) + "*" + str(
                        fragment["fragLength"])
                else:
                    nomeFilefine += str(re[1]) + str(re[2])
            if y != (len(lista) - 1) and i == 0:
                nomeFilefine += "_"
            elif y == (len(lista) - 1) and i == 0:
                if nModel != "":
                    nomeFilefine += "_" + nModel

            (consumedResidue, myHea) = actionImplemented(fragment["model"], fragment["chain"], re, consumedResidue,
                                                         myHea)
            # print "uuuuuuuuuuu",fragment["model"],fragment["chain"],re

    # print "fine processo residui delle strutture"

    # for uio in consumedResidue:
    #       print uio.get_id(), uio.get_resname()

    atom_number = 1
    previousChain = ""
    previousResName = ""
    previousResSeq = None
    lastRes = None
    previousIcode = ""
    pdbString = ""
    pdbString += "REMARK TITLE " + nomeFilefine + "\n"
    resiNumbering = 1
    dizioConvRes = {}
    for values, ite in myHea:
        item = ite[0]
        element = ite[1]
        orig_atom_num = item.get_serial_number()
        hetfield, resseq, icode = item.get_parent().get_id()
        segid = item.get_parent().get_segid()
        resname = item.get_parent().get_resname()
        chain_id = item.get_parent().get_parent().get_id()
        if previousResSeq != None:
            resaN = item.get_parent()
            prevResC = lastRes.get_parent()
            if not checkContinuity(prevResC, resaN) and resseq != previousResSeq:
                format_string = "%s%6i %-4s%c%3s %c%4i \n"
                if useDizioConv:
                    arg = ("TER  ", atom_number, ' ', ' ', previousResName, previousChain, resiNumbering)
                    pdbString += (format_string % arg)
                if len(dizioConv.keys()) == 0:
                    resiNumbering += 20  # to break the continuity of the residues
                    atom_number += 20
                    dizioConvRes[(chain_id, (hetfield, resseq, icode), element)] = (atom_number, resiNumbering)

        if previousResSeq != None and resseq != previousResSeq or (resseq == previousResSeq and icode != previousIcode):
            resiNumbering += 1

        if len(dizioConv.keys()) > 0:
            resiNumbering = (dizioConv[(chain_id, (hetfield, resseq, icode), element)])[1]
            atom_number = (dizioConv[(chain_id, (hetfield, resseq, icode), element)])[0]

        if not useDizioConv:
            pdbString += get_atom_line(item, element, hetfield, segid, orig_atom_num, resname, resseq, icode, chain_id)
        else:
            pdbString += get_atom_line(item, element, chr(ord(' ')), segid, atom_number, resname, resiNumbering,
                                     chr(ord(' ')), chain_id, normalize=normalize, bfactorNor=bfactorNor)
        if len(dizioConv.keys()) == 0:
            dizioConvRes[(chain_id, (hetfield, resseq, icode), element)] = (atom_number, resiNumbering)

        atom_number += 1
        previousChain = chain_id
        previousResName = resname
        previousResSeq = resseq
        previousIcode = icode
        lastRes = item
    if previousResSeq != None:
        format_string = "%s%6i %-4s%c%3s %c%4i \n"
        if useDizioConv:
            arg = ("TER  ", atom_number, ' ', ' ', previousResName, previousChain, resiNumbering)
            pdbString += (format_string % arg)
            atom_number += 1
    pdbString += "END\n"

    if len(dizioConv.keys()) == 0:
        return (nomeFilefine, pdbString, dizioConvRes)
    else:
        return (nomeFilefine, pdbString)

#TODO: This function is needed by getSuperimp but must be deleted with it
class PriorityEntry(object):

    def __init__(self, priority, data):
        self.priority = priority
        self.data = data

    def formatted(self):
        return [self.priority, self.data]

    def __lt__(self, other):
        return self.priority < other.priority

#TODO: This function is needed by getSuperimp but must be deleted with it
class Heap(object):
    """
    :author: Dr. Massimo Domenico Sammito
    :email: msacri@ibmb.csic.es / massimo.sammito@gmail.com

    """

    def __init__(self):
        """ create a new min-heap. """
        self._heap = []

    def push(self, priority, item):
        """ Push an item with priority into the heap.
            Priority 0 is the highest, which means that such an item will
            be popped first."""
        # assert priority >= 0
        aux = PriorityEntry(priority, item)
        heapq.heappush(self._heap, aux)

    def pop(self):
        """ Returns the item with lowest priority. """
        element = heapq.heappop(self._heap)  # (prio, item)[1] == item
        (prio, item) = element.formatted()
        return (prio, item)

    def len(self):
        return len(self._heap)

    def asList(self):
        aux = [x.formatted() for x in self._heap]
        return aux

    def __iter__(self):
        """ Get all elements ordered by asc. priority. """
        return self

    def __next__(self):
        """ Get all elements ordered by their priority (lowest first). """
        try:
            return self.pop()
        except IndexError:
            raise StopIteration

def change_chain(pdb, chain, atom_list=["ATOM  ", "ANISOU", "HETATM", "TER   "]):
    allpdb = pdb.splitlines()
    out = ""
    for line in allpdb:
        # only look at records indicated by atom_list
        if line[0:6] not in atom_list:
            if not line.startswith("END"):
                out += line + "\n"
            continue

        # Grab only residues belonging to chain
        out += line[:21] + chain + line[22:] + "\n"
    return out

#TODO: Used in compareDistributionAccordingOrientation
def getListValidFromTwoList(refe, compare, areFragments=False):
    listCombi = itertools.permutations(compare)
    listValid = []

    # print "reference",len(refe),"compare",len(compare)
    # printSecondaryStructureElements(refe)
    lio = []
    # print "--",len(compare)
    for combi in listCombi:
        lio.append(combi)
    listCombi = lio

    for combi in listCombi:
        valid = False
        # print combi
        for ui in range(len(refe)):
            # print refe
            frag1 = refe[ui]
            frag2 = combi[ui]
            valuta = None
            if areFragments:
                valuta = frag2["fragLength"] != frag1["fragLength"] and (
                    (frag2["sstype"] in ["bs", "cbs"] and frag1["sstype"] not in ["bs", "cbs"]) or (
                        frag2["sstype"] in ["ah", "ch"] and frag1["sstype"] not in ["ah", "ch"]))
            else:
                valuta = len(frag2) != len(frag1)

            if valuta:
                valid = False
                break
            else:
                valid = True
        if valid:
            listValid.append(combi)
    return listValid

#TODO: Used in SELSLIB2
def compareDistributionAccordingOrientation(AlistFrags, BlistFrags, threshold, shift, where, returnCVS=False):
    returnList = []
    nTota = []
    returnAnyway = []
    retcvs = []
    BlistFragsVali = getListValidFromTwoList(AlistFrags, BlistFrags, areFragments=True)
    bestDegree = 10000
    bestTuple = None

    for BlistFrags in BlistFragsVali:
        for i in range(len(AlistFrags)):
            refrag = []
            reAny = []
            recvs = []
            frag1 = AlistFrags[i]
            frag2 = BlistFrags[i]
            cvv1 = generateVectorsCVS(frag1, full=True)
            cvv2 = generateVectorsCVS(frag2, full=True)

            # print "---cvv1",len(cvv1)
            # print "---cvv2",len(cvv2)
            if len(cvv1) != len(cvv2):
                cvv1 = [cvv1[0], cvv1[int(len(cvv1) / 2)], cvv1[-1]]
                cvv2 = [cvv2[0], cvv2[int(len(cvv2) / 2)], cvv2[-1]]
                # continue

            nThis = len(cvv1)
            ind1, cv1, ind2, cv2, cav1, ov1, cav2, ov2 = 0, 0, 0, 0, 0, 0, 0, 0

            for t in range(len(cvv1) - shift):
                if shift != 0 and where == "A":
                    ind1, cv1, cav1, ov1 = cvv1[t + shift]
                    ind2, cv2, cav2, ov2 = cvv2[t]
                elif shift != 0 and where == "B":
                    ind1, cv1, cav1, ov1 = cvv1[t]
                    ind2, cv2, cav2, ov2 = cvv2[t + shift]
                elif shift == 0:
                    ind1, cv1, cav1, ov1 = cvv1[t]
                    ind2, cv2, cav2, ov2 = cvv2[t]

                X1 = cv1[0]
                Y1 = cv1[1]
                Z1 = cv1[2]
                X2 = cv2[0]
                Y2 = cv2[1]
                Z2 = cv2[2]
                TetaReal = angle_between([X1, Y1, Z1], [X2, Y2, Z2], [1.0, 1.0, 1.0], signed=False)
                tetadeg = TetaReal * 57.2957795

                if tetadeg <= threshold:
                    refrag.append((ind1, ind2, tetadeg))
                reAny.append((ind1, ind2, tetadeg))
                recvs.append((cav1, ov1, cav2, ov2))

            returnList.append(refrag)
            returnAnyway.append(reAny)
            nTota.append(nThis)
            retcvs.append(recvs)

        media = []
        for i in range(len(nTota)):
            summa = 0
            for t in range(len(returnAnyway[i])):
                summa += float(((returnAnyway[i])[t])[2])
            media.append(summa / len(returnAnyway[i]))
        vas = max(media)
        if vas < bestDegree:
            bestDegree = vas
            if not returnCVS:
                bestTuple = nTota, returnList, returnAnyway
            else:
                bestTuple = nTota, returnList, returnAnyway, retcvs

    return bestTuple

def generateVectorsCVS(fragment, full=False):
    distcvs = []
    nS = fragment["fragLength"] - SUFRAGLENGTH + 1
    CAatomCoord = fragment["CAatomCoord"]
    OatomCoord = fragment["OatomCoord"]

    for plo in range(nS):
        # print "plo is",plo
        xca = 0.0
        yca = 0.0
        zca = 0.0
        xo = 0.0
        yo = 0.0
        zo = 0.0
        for qlo in range(SUFRAGLENGTH):
            # print "\tqlo is",qlo
            xca += CAatomCoord[plo + qlo][0]
            yca += CAatomCoord[plo + qlo][1]
            zca += CAatomCoord[plo + qlo][2]
            xo += OatomCoord[plo + qlo][0]
            yo += OatomCoord[plo + qlo][1]
            zo += OatomCoord[plo + qlo][2]
        xca /= SUFRAGLENGTH
        yca /= SUFRAGLENGTH
        zca /= SUFRAGLENGTH
        xo /= SUFRAGLENGTH
        yo /= SUFRAGLENGTH
        zo /= SUFRAGLENGTH

        XH = xca - xo
        YH = yca - yo
        ZH = zca - zo
        if not full:
            distcvs.append((plo, [XH, YH, ZH]))
        else:
            distcvs.append((plo, [XH, YH, ZH], [xca, yca, zca], [xo, yo, zo]))

    return distcvs

def reading_HHPRED_file(hhr_file):

    """Given a multiple aligment file from HHPRED, parse it and generates a dictionary containing the information

    :author: Ana Medina
    :email: ambcri@ibmb.csic.es

    :param hhr_file: path to the .hhr file containing multiple alignment
    :type hhr_file: str
    :return: dict_homologs
    :rtype dict_homologs: dictionary
    """

    f = open(hhr_file, 'r')
    HHpredfile = f.readlines()
    f.close()

    homolog = ''
    new_pdb = True
    dic_homologs = {}
    for i, line in enumerate(HHpredfile):
        if line.startswith('>')and line[1:7] not in dic_homologs:
            homolog = line[1:7]
            dic_homologs[homolog] = {}
            dic_homologs[homolog]['sequence'] = []
            dic_homologs[homolog]['alignment'] = []
            dic_homologs[homolog]['start'] = []
            try:
                dic_homologs[homolog]['related'] = line.split('Related PDB entries: ')[1].split()
            except:
                dic_homologs[homolog]['related'] = []
            new_pdb = False
            seq_homolog = ''
        elif line.startswith('>')and line[1:7] in dic_homologs:
            new_pdb = True
            homolog = ''
        elif homolog and not new_pdb and line.startswith('T Consensus'):
            sequence = ''
            align = ''
            for e,letter in enumerate(HHpredfile[i + 1].strip()[22:-11]):
                if letter != '-' and letter != ' ':
                    sequence += letter
                    #print(HHpredfile[i - 1][22:], letter, e, HHpredfile[i - 1][22:][e])
                    align += (HHpredfile[i - 1][22:])[e]
            dic_homologs[homolog]['sequence'].append(sequence)
            dic_homologs[homolog]['alignment'].append(align)
            dic_homologs[homolog]['start'].append(int(HHpredfile[i].split()[2]))
    return (dic_homologs)


def get_sequence_from_pdb(pdb):
    """Given a pdb produces a dictinary where each key is a chain. The value is another dictionary that contains
    two different keys: sequence in one letter code and starting residue number of this chain.
    If missing aminoacids according to residue numbering, X is set in sequence.

    :author: Ana Medina
    :email: ambcri@ibmb.csic.es

    :param pdb: paths to the pdb to read
    :type pdb: str
    :return: dict_sequence: dictionary with the sequence for each chain
    :rtype dict_sequence: dictionary
    """

    f = open (pdb, 'r')
    lines = f.readlines()
    f.close()

    dict_sequence = {}
    sequence = ''

    atoms = [line for line in lines if line.startswith('ATOM')]

    for i, line in enumerate(atoms):
        current_res = int(line[22:26].strip())
        res = line[17:20]
        current_chain = line[20:22].strip()
        if i == 0:
            last_chain = line[20:22].strip()
            first_res = int(line[22:26].strip())
            dict_sequence[last_chain] = {}

        elif current_res == last_res or res not in AADICMAP:
            continue

        elif current_chain != last_chain:
            if not current_chain in dict_sequence:
                dict_sequence[last_chain]['sequence'] = sequence
                dict_sequence[last_chain]['first_res'] = first_res
                last_chain = current_chain
                dict_sequence[last_chain] = {}
                first_res = int(line[22:26].strip())
                sequence = ''

            else:
                sys.exit()

        elif current_res > last_res +1:
            sequence+=('X' * (current_res - last_res - 1))

        last_res = current_res
        try:
            sequence+=(AADICMAP[res])
        except:
            continue

    dict_sequence[last_chain]['sequence'] = sequence
    dict_sequence[last_chain]['first_res'] = first_res

    return dict_sequence

def prepare_model_for_spack (model_annotated, crystcard):
    f = open(model_annotated, 'r')
    lines = f.readlines()
    f.close()

    input_spack = os.path.join(model_annotated[:-4] + '_spack.pdb')
    g = open(input_spack, 'w') #Removing .pdb termination and adding _spack.pdb
    g.write(crystcard)

    for line in lines:
        if line.startswith('ATOM'):
            g.write(line[:16] + ' ' + line[17:56] + '1.00' + line[60:]) #In case of disordered atoms
        else:
            g.write(line)
    g.close()

    return (input_spack)

def get_crystcard(reference):
    crystcard = None

    f = open(reference, "r")
    all_lines = f.readlines()
    f.close()


    for line in all_lines:
        if line[:6] in ["CRYST1"]:
            crystcard = line
            # pdb_search_in += line
        # if line[:6] in ["TITLE ", "SCALE1", "SCALE2", "SCALE3"]:
        #     pdb_search_in += line

    if not crystcard:
        raise Exception("ERROR: There was not CRYST1 remark present in the input model. Symmetry copies cannot be generated")

    return crystcard

def extract_cryst_card_pdb(pdb):
    """ Extract a cryst1 from a pdb file
    :param pdb: PDB
    :type pdb: str (Can be a file or a string with the pdb information)
    :rtype: None or crystalcard
    """
    if os.path.isfile(pdb):
        with open(pdb, 'r') as handle:
            pdb = handle.read()

    pdb_lines = pdb.split("\n")
    for linea in pdb_lines:
        if linea.startswith("CRYST1"):  # E.G. CRYST1   30.279   91.989   32.864  90.00 112.60  90.00 P 1 21 1      2
            cryst_card = linea
            return cryst_card
    return None

def add_cryst_card_pdb(pdb_file, cryst_card=None, space_group=None, cell_dim=None):
    """ Add a cryst1 record to a pdb file
    :param path_pdb: path to the pdb file to modify. Note that it will be overwritten
    :type path_pdb: str
    :param cryst_card: cryst1 record
    :type cryst_card: str
    :param space_group: space_group
    :type space_group: str
    :param cell_dim: cell_dim record
    :type cell_dim: str
    :return: True for success, False otherwise
    :rtype: bool
    """
    try:
        if cryst_card is None:
            CRD = '{r:<6}{a:>9.3f}{b:>9.3f}{c:>9.3f}{alpha:7.2f}{beta:7.2f}{gamma:7.2f} {sp:<11}{z:>4}'
            data = dict(r='CRYST1', a=float(cell_dim[0]), b=float(cell_dim[1]), c=float(cell_dim[2]), alpha=float(cell_dim[3]),
                        beta=float(cell_dim[4]), gamma=float(cell_dim[5]), sp=space_group, z=1)
            cryst_card = CRD.format(**data)
        with open(pdb_file, 'r') as handle:
            pdb_dump = handle.read()
        with open(pdb_file, 'w') as handle:
            handle.write(cryst_card + "\n")
            handle.write(pdb_dump)
        return True
    except:
        print('ERROR: Something went wrong adding the CRYST1 record to the pdb at', pdb_file)
        return False

def extract_fasta_pdb(pdb_file, NC=1, fasta_file=None):
    """ Get the fasta from a pdb file
    :param pdb_file: PDB
    :type pdb_file: str
    :param fasta_file: Path to a file
    :type fasta_file: str
    :rtype: None or FASTA
    """
    try:
        record = SeqIO.read(pdb_file, "pdb-atom")
        if fasta_file is not None:
            with open(fasta_file, "w+") as handle:
                for i in range(NC):
                    SeqIO.write(record, handle, "fasta")
        return record.seq
    except:
        print('ERROR: Something went wrong extracting the fasta record from the pdb at', pdb_file)
        pass
    return None




