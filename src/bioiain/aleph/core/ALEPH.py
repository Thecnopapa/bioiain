#! /usr/bin/env python
# -*- coding: utf-8 -*-

#future imports
from __future__ import print_function
from __future__ import division
from future.standard_library import install_aliases
install_aliases()
#from __future__ import unicode_literals

# System imports
import os
import sys
import stat
import typing

#####################################################################
# All the functions in this module were written by:                 #
# author: Massimo Sammito                                           #
# email: massimo.sammito@gmail.com                                  #
#####################################################################


info_p = sys.version_info
info_g = (sys.version).splitlines()
PYTHON_V = info_p.major

import argparse
import warnings
warnings.filterwarnings("ignore")

import time
import datetime
import copy
import traceback
import shutil
import subprocess
import multiprocessing
import io
import tempfile
from io import BytesIO
import xml.etree.ElementTree as ET

#System imports
from operator import itemgetter

if PYTHON_V == 3:
    import configparser
    from urllib.request import urlopen, Request
    from builtins import range
    from builtins import str


# format imports
import json
import pickle
import tarfile
import gzip
import shelve

# graphics
import matplotlib

matplotlib.use('Agg')
from termcolor import colored, cprint
import matplotlib.pyplot as plt
import matplotlib.transforms as transforms
from matplotlib.lines import Line2D

# Scientific and numerical imports
import math
import numpy
import random
import scipy
import scipy.spatial
import scipy.sparse
import scipy.cluster
import scipy.special
import igraph
from igraph import *
#import igraph.vendor.texttable
import Bio.PDB
from Bio import pairwise2
from Bio.Align import substitution_matrices
names = substitution_matrices.load()
#from Bio import SeqIO
#import csb.bio.utils

import itertools
from itertools import permutations
from itertools import combinations

import pandas
from sklearn.cluster import AgglomerativeClustering

# personal libraries imports
from . import Bioinformatics
from . import SystemUtility
from . import Grid

warnings.simplefilter("ignore", Bio.PDB.PDBExceptions.PDBConstructionWarning)
#warnings.simplefilter("ignore", DeprecationWarning)

#######################################################################################################
#                                           CONSTANT VARIABLES                                         #
#######################################################################################################

color_dict = {}
list_ids = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz1234567890"
kn = list(igraph.drawing.colors.known_colors.keys())
numpy.random.shuffle(kn)
for h, idn in enumerate(list_ids):
    color_dict[idn] = kn[h]

SCALING = "min_max"
# :(0.5,0.1,0.3),"min_max":(0.5,0.1,0.3)}
THRESH_CORR = {"robust_scale":{"remove_redundance_chain":0.0005,"compare_instruction":0.005},
               "min_max":     {"remove_redundance_chain":0.0005,"compare_instruction":0.005}}

BS_UU_EA = {0: 24.6315759411, 1: 24.6315759411, 2: 24.6315759411, 3: 24.6315759411, 4: 24.6315759411, 5: 24.6315759411, 6: 24.6315759411, 7: 24.6315759411, 8: 24.6315759411, 9: 24.6315759411, 10: 25.181356029, 11: 25.181356029, 12: 25.181356029, 13: 25.181356029, 14: 25.181356029, 15: 25.181356029, 16: 25.181356029, 17: 25.181356029, 18: 25.181356029, 19: 25.181356029, 20: 19.4507911121, 21: 19.4507911121, 22: 19.4507911121, 23: 19.4507911121, 24: 19.4507911121, 25: 19.4507911121, 26: 19.4507911121, 27: 19.4507911121, 28: 19.4507911121, 29: 19.4507911121, 30: 11.3597418176, 31: 11.3597418176, 32: 11.3597418176, 33: 11.3597418176, 34: 11.3597418176, 35: 11.3597418176, 36: 11.3597418176, 37: 11.3597418176, 38: 11.3597418176, 39: 11.3597418176, 40: 7.02861712458, 41: 7.02861712458, 42: 7.02861712458, 43: 7.02861712458, 44: 7.02861712458, 45: 7.02861712458, 46: 7.02861712458, 47: 7.02861712458, 48: 7.02861712458, 49: 7.02861712458, 50: 4.80950476952, 51: 4.80950476952, 52: 4.80950476952, 53: 4.80950476952, 54: 4.80950476952, 55: 4.80950476952, 56: 4.80950476952, 57: 4.80950476952, 58: 4.80950476952, 59: 4.80950476952, 60: 2.9131204661, 61: 2.9131204661, 62: 2.9131204661, 63: 2.9131204661, 64: 2.9131204661, 65: 2.9131204661, 66: 2.9131204661, 67: 2.9131204661, 68: 2.9131204661, 69: 2.9131204661, 70: 2.2662363626, 71: 2.2662363626, 72: 2.2662363626, 73: 2.2662363626, 74: 2.2662363626, 75: 2.2662363626, 76: 2.2662363626, 77: 2.2662363626, 78: 2.2662363626, 79: 2.2662363626, 80: 1.00816816131, 81: 1.00816816131, 82: 1.00816816131, 83: 1.00816816131, 84: 1.00816816131, 85: 1.00816816131, 86: 1.00816816131, 87: 1.00816816131, 88: 1.00816816131, 89: 1.00816816131, 90: 0.658308105329, 91: 0.658308105329, 92: 0.658308105329, 93: 0.658308105329, 94: 0.658308105329, 95: 0.658308105329, 96: 0.658308105329, 97: 0.658308105329, 98: 0.658308105329, 99: 0.658308105329, 100: 0.26275204204, 101: 0.26275204204, 102: 0.26275204204, 103: 0.26275204204, 104: 0.26275204204, 105: 0.26275204204, 106: 0.26275204204, 107: 0.26275204204, 108: 0.26275204204, 109: 0.26275204204, 110: 0.145656023305, 111: 0.145656023305, 112: 0.145656023305, 113: 0.145656023305, 114: 0.145656023305, 115: 0.145656023305, 116: 0.145656023305, 117: 0.145656023305, 118: 0.145656023305, 119: 0.145656023305, 120: 0.129948020792, 121: 0.129948020792, 122: 0.129948020792, 123: 0.129948020792, 124: 0.129948020792, 125: 0.129948020792, 126: 0.129948020792, 127: 0.129948020792, 128: 0.129948020792, 129: 0.129948020792, 130: 0.0842520134803, 131: 0.0842520134803, 132: 0.0842520134803, 133: 0.0842520134803, 134: 0.0842520134803, 135: 0.0842520134803, 136: 0.0842520134803, 137: 0.0842520134803, 138: 0.0842520134803, 139: 0.0842520134803, 140: 0.0514080082253, 141: 0.0514080082253, 142: 0.0514080082253, 143: 0.0514080082253, 144: 0.0514080082253, 145: 0.0514080082253, 146: 0.0514080082253, 147: 0.0514080082253, 148: 0.0514080082253, 149: 0.0514080082253, 150: 0.0114240018278, 151: 0.0114240018278, 152: 0.0114240018278, 153: 0.0114240018278, 154: 0.0114240018278, 155: 0.0114240018278, 156: 0.0114240018278, 157: 0.0114240018278, 158: 0.0114240018278, 159: 0.0114240018278, 160: 0.0071400011424, 161: 0.0071400011424, 162: 0.0071400011424, 163: 0.0071400011424, 164: 0.0071400011424, 165: 0.0071400011424, 166: 0.0071400011424, 167: 0.0071400011424, 168: 0.0071400011424, 169: 0.0071400011424, 170: 0.0, 171: 0.0, 172: 0.0, 173: 0.0, 174: 0.0, 175: 0.0, 176: 0.0, 177: 0.0, 178: 0.0, 179: 0.0, 180: 0.0}
BS_UD_EA = {0: 0.00200625952973, 1: 0.00200625952973, 2: 0.00200625952973, 3: 0.00200625952973, 4: 0.00200625952973, 5: 0.00200625952973, 6: 0.00200625952973, 7: 0.00200625952973, 8: 0.00200625952973, 9: 0.00200625952973, 10: 0.0170532060027, 11: 0.0170532060027, 12: 0.0170532060027, 13: 0.0170532060027, 14: 0.0170532060027, 15: 0.0170532060027, 16: 0.0170532060027, 17: 0.0170532060027, 18: 0.0170532060027, 19: 0.0170532060027, 20: 0.00702190835406, 21: 0.00702190835406, 22: 0.00702190835406, 23: 0.00702190835406, 24: 0.00702190835406, 25: 0.00702190835406, 26: 0.00702190835406, 27: 0.00702190835406, 28: 0.00702190835406, 29: 0.00702190835406, 30: 0.0060187785892, 31: 0.0060187785892, 32: 0.0060187785892, 33: 0.0060187785892, 34: 0.0060187785892, 35: 0.0060187785892, 36: 0.0060187785892, 37: 0.0060187785892, 38: 0.0060187785892, 39: 0.0060187785892, 40: 0.00702190835406, 41: 0.00702190835406, 42: 0.00702190835406, 43: 0.00702190835406, 44: 0.00702190835406, 45: 0.00702190835406, 46: 0.00702190835406, 47: 0.00702190835406, 48: 0.00702190835406, 49: 0.00702190835406, 50: 0.0230719845919, 51: 0.0230719845919, 52: 0.0230719845919, 53: 0.0230719845919, 54: 0.0230719845919, 55: 0.0230719845919, 56: 0.0230719845919, 57: 0.0230719845919, 58: 0.0230719845919, 59: 0.0230719845919, 60: 0.0290907631811, 61: 0.0290907631811, 62: 0.0290907631811, 63: 0.0290907631811, 64: 0.0290907631811, 65: 0.0290907631811, 66: 0.0290907631811, 67: 0.0290907631811, 68: 0.0290907631811, 69: 0.0290907631811, 70: 0.0642003049514, 71: 0.0642003049514, 72: 0.0642003049514, 73: 0.0642003049514, 74: 0.0642003049514, 75: 0.0642003049514, 76: 0.0642003049514, 77: 0.0642003049514, 78: 0.0642003049514, 79: 0.0642003049514, 80: 0.155485113554, 81: 0.155485113554, 82: 0.155485113554, 83: 0.155485113554, 84: 0.155485113554, 85: 0.155485113554, 86: 0.155485113554, 87: 0.155485113554, 88: 0.155485113554, 89: 0.155485113554, 90: 0.489527325255, 91: 0.489527325255, 92: 0.489527325255, 93: 0.489527325255, 94: 0.489527325255, 95: 0.489527325255, 96: 0.489527325255, 97: 0.489527325255, 98: 0.489527325255, 99: 0.489527325255, 100: 1.53880105931, 101: 1.53880105931, 102: 1.53880105931, 103: 1.53880105931, 104: 1.53880105931, 105: 1.53880105931, 106: 1.53880105931, 107: 1.53880105931, 108: 1.53880105931, 109: 1.53880105931, 110: 4.22117005056, 111: 4.22117005056, 112: 4.22117005056, 113: 4.22117005056, 114: 4.22117005056, 115: 4.22117005056, 116: 4.22117005056, 117: 4.22117005056, 118: 4.22117005056, 119: 4.22117005056, 120: 10.0764384881, 121: 10.0764384881, 122: 10.0764384881, 123: 10.0764384881, 124: 10.0764384881, 125: 10.0764384881, 126: 10.0764384881, 127: 10.0764384881, 128: 10.0764384881, 129: 10.0764384881, 130: 14.2544739588, 131: 14.2544739588, 132: 14.2544739588, 133: 14.2544739588, 134: 14.2544739588, 135: 14.2544739588, 136: 14.2544739588, 137: 14.2544739588, 138: 14.2544739588, 139: 14.2544739588, 140: 18.0824171415, 141: 18.0824171415, 142: 18.0824171415, 143: 18.0824171415, 144: 18.0824171415, 145: 18.0824171415, 146: 18.0824171415, 147: 18.0824171415, 148: 18.0824171415, 149: 18.0824171415, 150: 17.8075595859, 151: 17.8075595859, 152: 17.8075595859, 153: 17.8075595859, 154: 17.8075595859, 155: 17.8075595859, 156: 17.8075595859, 157: 17.8075595859, 158: 17.8075595859, 159: 17.8075595859, 160: 18.9531337774, 161: 18.9531337774, 162: 18.9531337774, 163: 18.9531337774, 164: 18.9531337774, 165: 18.9531337774, 166: 18.9531337774, 167: 18.9531337774, 168: 18.9531337774, 169: 18.9531337774, 170: 14.2655083862, 171: 14.2655083862, 172: 14.2655083862, 173: 14.2655083862, 174: 14.2655083862, 175: 14.2655083862, 176: 14.2655083862, 177: 14.2655083862, 178: 14.2655083862, 179: 14.2655083862, 180: 14.2655083862}
BS_MAX = {0:max(BS_UD_EA.values()), 1:max(BS_UU_EA.values())}

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

dict_similarity = {1: 45, 2: 56, 3: 67, 4: 78, 5: 89}

PATH_NEW_BORGES = ""
PATH_PYTHON_INTERPRETER = ""
GRID_TYPE = ""
MAX_PDB_TAR = 5.0
NUMBER_OF_PARALLEL_GRID_JOBS = 70000
MAX_EVALUE = 0.005
number_of_solutions = 0
canCluster = True
doCluster_global = True
listjobs = []
toclientdir = ""

#######################################################################################################
#                                           SUPPORT FUNC                                              #
#######################################################################################################

def get_blast_models(sequence, max_evalue=0.005):
    search_string1 = "https://search.rcsb.org/rcsbsearch/v1/query?json=%7B%22query%22:%7B%22type%22:%22terminal%22,%22service%22:%22sequence%22,%22parameters%22:%7B%22evalue_cutoff%22:+1,%22target%22:%22pdb_protein_sequence%22,%22value%22:%22" + sequence.upper() + "%22%7D%7D,%22request_options%22:%7B%22scoring_strategy%22:+%22sequence%22%7D,%22return_type%22:+%22polymer_entity%22%7D"
    search_string2 = "https://search.rcsb.org/rcsbsearch/v1/query?json=%7B%22query%22:%7B%22type%22:%22terminal%22,%22service%22:%22sequence%22,%22parameters%22:%7B%22evalue_cutoff%22:+1,%22target%22:%22pdb_protein_sequence%22,%22value%22:%22" + sequence.upper() + "%22%7D%7D,%22request_options%22:%7B%22scoring_strategy%22:+%22sequence%22%7D,%22return_type%22:+%22polymer_instance%22%7D"

    ebiresp = SystemUtility.request_url(search_string2, "www.rcsb.org")
    if not ebiresp: return {}
    infol = json.loads(ebiresp)
    dictio = {}
    for q, node in enumerate(infol["result_set"]):
        pdbid = node["identifier"].split(".")[0].lower()
        if pdbid not in dictio:
            dictio[pdbid] = [node["identifier"].split(".")[1]]
        else:
            dictio[pdbid].append(node["identifier"].split(".")[1])

    ebiresp = SystemUtility.request_url(search_string1, "www.rcsb.org")
    if not ebiresp: return {}
    infol = json.loads(ebiresp)
    hits = {q + 1: {"assemblies": {},
                    "pdbid": node["identifier"].split("_")[0].lower(),
                    "chains": dictio[node["identifier"].split("_")[0].lower()],
                    "score": float(node["score"]),
                    "evalue": float(node["services"][0]["nodes"][0]["match_context"][0]["evalue"]),
                    "identity": float(node["services"][0]["nodes"][0]["match_context"][0]["sequence_identity"]),
                    "lenalign": int(node["services"][0]["nodes"][0]["match_context"][0]["alignment_length"]),
                    "target_size": int(node["services"][0]["nodes"][0]["match_context"][0]["subject_length"])
                    } for q, node in enumerate(infol["result_set"]) if
            node["identifier"].split("_")[0].lower() in dictio}

    hits = {key: value for key, value in hits.items() if value["evalue"] <= max_evalue}
    # print("=========================HITS==============================")
    # print(hits)
    # print("===========================================================")
    return hits

def get_blast_models_old(sequence):
    global MAX_EVALUE

    search_string = "https://www.rcsb.org/pdb/rest/getBlastPDB1?sequence=%s&eCutOff=10.0&matrix=BLOSUM62&outputFormat=XML" % (sequence.upper())
    req = Request(search_string)
    z = urlopen(req).read()
    f = BytesIO(z)
    tree = ET.parse(f)
    root = tree.getroot()
    nodes = root.findall('BlastOutput_iterations/Iteration/Iteration_hits/Hit')
    hits = {int(node.find("Hit_num").text): {"assemblies": {},
                                             "pdbid": node.find("Hit_def").text.split("|")[0].split(":")[0].lower(),
                                             "chains": node.find("Hit_def").text.split("|")[0].split(":")[2].split(","),
                                             "evalue": float(node.find("Hit_hsps/Hsp/Hsp_evalue").text),
                                             "identity": int(node.find("Hit_hsps/Hsp/Hsp_identity").text),
                                             "lenalign": int(node.find("Hit_hsps/Hsp/Hsp_align-len").text),
                                             "target_size": int(node.find("Hit_len").text)} for node in nodes}
    hits = {key: value for key, value in hits.items() if value["evalue"] <= MAX_EVALUE}

    return hits

def get_structure_ids_from_pdbids(pdbids):
    pfams = []
    cathids = []
    scops = []
    for value in pdbids:
        pdbid,chains = value
        search_string = "https://www.ebi.ac.uk/pdbe/api/mappings/%s" % (pdbid.lower())
        req = Request(search_string)
        z = urlopen(req).read()
        dicr = json.loads(z)
        pfams += [key for key in dicr[pdbid.lower()]["Pfam"].keys() if len(set([oo["chain_id"] for oo in dicr[pdbid.lower()]["Pfam"][key]["mappings"]])&set(chains))>0]
        cathids += [key for key in dicr[pdbid.lower()]["CATH"].keys() if len(set([oo["chain_id"] for oo in dicr[pdbid.lower()]["CATH"][key]["mappings"]])&set(chains))>0]
        scops += [dicr[pdbid.lower()]["SCOP"][key]["superfamily"]["description"] for key in dicr[pdbid.lower()]["SCOP"].keys() if "superfamily" in dicr[pdbid.lower()]["SCOP"][key] and len(set([oo["chain_id"] for oo in dicr[pdbid.lower()]["SCOP"][key]["mappings"]])&set(chains))>0]
    return list(set(pfams)),list(set(cathids)),list(set(scops))

def download_database(pfams=None,cathids=None,scops=None):
    if not pfams and not cathids and not scops: return None
    pdbs = []
    if pfams is not None:
        query = "("+"%20OR%20".join([pfam for pfam in pfams])+")"
        search_string = "https://www.ebi.ac.uk/pdbe/search/pdb/select?q=pfam_accession:%s&wt=json&fl=pdb_id,homologus_pdb_entity_id&rows=10000" % (query)
        #print(search_string)
        req = Request(search_string)
        z = urlopen(req).read()
        dicr = json.loads(z)
        pdbs += [doc["pdb_id"].lower() for doc in dicr["response"]["docs"] if "pdb_id" in doc]+[t for x in [[pdbr.split("_")[0] for pdbr in doc["homologus_pdb_entity_id"]] for doc in dicr["response"]["docs"] if  "homologus_pdb_entity_id" in doc] for t in x]
    if cathids is not None:
        query = "(" + "%20OR%20".join([cath for cath in cathids]) + ")"
        search_string = "https://www.ebi.ac.uk/pdbe/search/pdb/select?q=cath_code:%s&wt=json&fl=pdb_id,homologus_pdb_entity_id&rows=10000" % (query)
        #print(search_string)
        req = Request(search_string)
        z = urlopen(req).read()
        dicr = json.loads(z)
        pdbs += [doc["pdb_id"].lower() for doc in dicr["response"]["docs"] if "pdb_id" in doc]+[t for x in [[pdbr.split("_")[0] for pdbr in doc["homologus_pdb_entity_id"]] for doc in dicr["response"]["docs"] if  "homologus_pdb_entity_id" in doc] for t in x]
    if scops is not None:
        query = "(" + "%20OR%20".join(["%20".join(scop.split()) for scop in scops]) + ")"
        search_string = "https://www.ebi.ac.uk/pdbe/search/pdb/select?q=scop_superfamily:%s&wt=json&fl=pdb_id,homologus_pdb_entity_id&rows=10000" % (
            query)
        print(search_string)
        req = Request(search_string)
        z = urlopen(req).read()
        dicr = json.loads(z)
        pdbs += [doc["pdb_id"].lower() for doc in dicr["response"]["docs"] if "pdb_id" in doc]+[t for x in [[pdbr.split("_")[0] for pdbr in doc["homologus_pdb_entity_id"]] for doc in dicr["response"]["docs"] if  "homologus_pdb_entity_id" in doc] for t in x]

    pdbs = list(set(pdbs))
    databse_from_seq = "./database_from_sequence"
    if not os.path.exists(databse_from_seq): os.makedirs(databse_from_seq)
    urlbase = "https://files.rcsb.org/download/"
    for pdb in pdbs:
        subdir = os.path.join(databse_from_seq,pdb[1:3])
        if not os.path.exists(subdir): os.makedirs(subdir)
        pdbf = os.path.join(subdir,str(pdb)+".pdb")
        if not os.path.exists(pdbf):
            try:
                req = Request(urlbase+str(pdb)+".pdb")
                z = urlopen(req).read()
                with open(pdbf,"wb") as f: f.write(z)
            except:
                print("The structure ",pdb,"is not available")

    return databse_from_seq

def consecutive_groups(iterable, ordering=lambda x: x):
    """Yield groups of consecutive items using :func:`itertools.groupby`.
    The *ordering* function determines whether two items are adjacent by
    returning their position.
    By default, the ordering function is the identity function. This is
    suitable for finding runs of numbers:
        >>> iterable = [1, 10, 11, 12, 20, 30, 31, 32, 33, 40]
        >>> for group in consecutive_groups(iterable):
        ...     print(list(group))
        [1]
        [10, 11, 12]
        [20]
        [30, 31, 32, 33]
        [40]
    For finding runs of adjacent letters, try using the :meth:`index` method
    of a string of letters:
        >>> from string import ascii_lowercase
        >>> iterable = 'abcdfgilmnop'
        >>> ordering = ascii_lowercase.index
        >>> for group in consecutive_groups(iterable, ordering):
        ...     print(list(group))
        ['a', 'b', 'c', 'd']
        ['f', 'g']
        ['i']
        ['l', 'm', 'n', 'o', 'p']
    """
    for k, g in itertools.groupby(
        enumerate(iterable), key=lambda x: x[0] - ordering(x[1])
    ):
        yield map(itemgetter(1), g)

class MinMaxScaler(object):
    def __init__(self, OldMin, OldMax, NewMin, NewMax):
        self.OldMin = OldMin
        self.OldMax = OldMax
        self.NewMin = NewMin
        self.NewMax = NewMax

    def scale(self, value, integer=False):
        q = (((value - self.OldMin) * (self.NewMax - self.NewMin)) / (self.OldMax - self.OldMin)) + self.NewMin
        if integer: return int(round(q))
        else: return q

#######################################################################################################
#                                           GEOMETRICAL FUNC                                          #
#######################################################################################################

def get_dot_between(A, B):
    """
    Return the dot product between two vectors
    :param A: A vector
    :type A: numpy.array
    :param B: B vector
    :type B: numpy.array
    :return dot: dot product
    :rtype: float
    """

    return (A * B).sum(axis=0)

def norm(A):
    return numpy.sqrt((A ** 2).sum(axis=0))

# @SystemUtility.profileit
# @SystemUtility.timing
def angle_between(A, B, N, signed=True):
    """
    ANGLE BETWEEN TWO 3D VECTORS:
    1- dot(norm(A),norm(B)) (ANGLES UNSIGNED, PROBLEMS FOR SMALL ANGLES WITH ROUNDINGS)
    2- arcos(dot(A,B)/(|A|*|B|))  (ANGLE UNSIGNED, PROBLEMS FOR SMALL ANGLES WITH ROUNDINGS)
    3- arctan2(|cross(A,B)|,dot(A,B)) (ANGLE UNSIGNED BUT NOT PROBLEMS OF ROUNDINGS
    Method 2 and 3 are equivalent: TESTED
      define a vector NORM ex.: N = [0,0,1]
      sign = dot(NORM,cross(A,B))
      if sign < 0 then ANGLE measured in 3 should be negative

    :param A: array coordinate first vector
    :type A: numpy.array
    :param B: array coordinate second vector
    :type B: numpy.array
    :param N: array coordinate normal direction
    :type N: numpy.array
    :param signed: True if signed angle should be returned, False otherwise
    :type signed: bool
    :return angle: the angle between A and B
    :rtype angle: float
    """

    # Cross = numpy.cross(A,B) ### this is 5 times slower according cProfiler,

    CrossX = A[1] * B[2] - A[2] * B[1]
    CrossY = A[2] * B[0] - A[0] * B[2]
    CrossZ = A[0] * B[1] - A[1] * B[0]
    Cross = numpy.asarray([CrossX, CrossY, CrossZ])

    fCross = numpy.sqrt((Cross ** 2).sum(axis=0))
    scaP2 = get_dot_between(A, B)
    Teta_2 = numpy.arctan2(fCross, scaP2)

    if signed:
        sign = get_dot_between(N, Cross)
        if sign < 0:
            Teta_2 = -Teta_2

        return Teta_2
    else:
        return Teta_2

def angle_between_by_acos(A, B):
    return numpy.arccos(get_dot_between(A, B) / (norm(A) * norm(B)))

def get_atoms_distance(in_a, in_b):
    """
    Return the euclidean distance between the two 3d coordinates
    :param in_a: 3d coord a
    :type in_a: numpy.array
    :param in_b: 3d coord b
    :type in_b: numpy.array
    :return d,D: distance and vector distance
    :rtype d,D: (float,numpy.array)
    """

    D = in_a - in_b
    d = numpy.sqrt((D ** 2).sum(axis=0))
    return d, D

#######################################################################################################
#                                          STARTER FUNCTIONS                                          #
#######################################################################################################

def perform_superposition_starter(args):
    global PYTHON_V
    start_time = time.time()

    Bioinformatics.rename_hetatm_and_icode(args.reference)
    if args.target and os.path.exists(args.target): Bioinformatics.rename_hetatm_and_icode(args.target)

    if os.path.exists(args.reference):
        dbtemp = os.path.join(tempfile.gettempdir(),"tempdirbase")
        if os.path.exists(dbtemp): shutil.rmtree(dbtemp)
        os.makedirs(dbtemp)
        shutil.copyfile(args.reference, os.path.join(dbtemp, os.path.basename(args.reference)))
        paths = []

        if args.targets and os.path.exists(args.targets):
            for root, subFolders, files in os.walk(args.targets):
                for fileu in files:
                    pdbf = os.path.join(root, fileu)
                    if pdbf.endswith(".pdb"):
                        paths.append(pdbf)
                        Bioinformatics.rename_hetatm_and_icode(pdbf)
        elif args.target and os.path.exists(args.target):
            paths.append(args.target)

        if not args.score_intra_fragment:
            try:
                args.score_intra_fragment = dict_similarity[args.similarity_intra_chain]
            except:
                args.score_intra_fragment = dict_similarity[3]
        if not args.score_inter_fragments:
            try:
                args.score_inter_fragments = dict_similarity[args.similarity_inter_chains]
            except:
                args.score_inter_fragments = dict_similarity[3]

        basename = os.path.basename(args.reference)[:-4]
        dic_models = {}

        dic_json = {'mode': 'superposition', 'submode': 'superposition based on Characteristic Vectors',
                    'parameterization': {'reference': args.reference,
                                         'strictness_ah': args.strictness_ah,
                                         'strictness_bs': args.strictness_bs, 'peptide_length': args.peptide_length,
                                         'score_intra_fragment': args.score_intra_fragment,
                                         'score_inter_fragments': args.score_inter_fragments,
                                         'criterium_selection_core': args.criterium_selection_core,
                                         'superposition rmsd threshold':args.rmsd_thresh, 'reverse':args.reverse}}
        modify_json_file('output.json', dic_json)

        for path in paths:
            generate_library(directory_database=dbtemp,
                             c_angle=-1, c_dist=-1, c_angle_dist=-1,
                             c_cvl_diff=-1,
                             score_intra_fragment=args.score_intra_fragment, j_angle=-1, j_dist=-1,
                             j_angle_dist=-1,
                             j_cvl_diff=-1, score_inter_fragments=args.score_inter_fragments,
                             rmsd_clustering=-1.0,
                             exclude_residues_superpose=1,
                             pdbmodel=path, weight="distance_avg",
                             strictness_ah=0.0001, strictness_bs=0.0001,
                             peptide_length=3, superpose=True, process_join=False,
                             nilges=args.nilges, ncssearch=True,
                             rmsd_max=args.rmsd_thresh, representative=args.representative,
                             swap_superposition=True, reverse=args.reverse, superposition_mode=True, multimer = args.multimer, interface = args.interface,
                             sidechains=args.sidechains, verbose=args.verbose)

            while len(multiprocessing.active_children()) > 0: pass

        while 1:
            if len(multiprocessing.active_children()) == 0:
                if args.reverse:
                    for root, subFolders, files in os.walk("./library"):
                        for fileu in files:
                            if fileu.startswith('output') and fileu.endswith('.pdb'):
                                pdbf = os.path.join(root, fileu)
                                shutil.move(pdbf, os.path.join(root, os.path.basename(args.reference)[:-4] + os.path.basename(pdbf)[6:]))

                if os.path.exists("./library/list_rmsd.txt"):
                    dic_models = {}
                    dic_json = {}
                    with open("./library/list_rmsd.txt", "r") as f0:
                        lines = f0.readlines()
                    with open("./library/list_rmsd.txt", "w") as f1:
                        for line in lines:
                            line_sep = line.split()
                            if args.reverse:
                                f1.write(os.path.basename(args.reference)[:-4] + line_sep[0][6:] + '\t' + line_sep[1] + '\t' + line_sep[2] + '\t' + line_sep[3] + '\n')
                                dic_models[os.path.basename(args.reference)[:-4] + line_sep[0][6:-4]] = [float(line_sep[1]), float(line_sep[2])]
                            else:
                                f1.write(line)
                                dic_models[line_sep[0][:-4]] = [float(line_sep[1]), float(line_sep[2])]

                        for key in dic_models:
                            dic_json[str(key) + '.pdb'] = {'rmsd and core': str(dic_models[key][0]) + "  " + str(dic_models[key][1])}
                            modify_json_file('output.json', dic_json, superposition=True)

                    dic_json = {"superposition": "{:.2f}s".format(time.time() - start_time)}
                    modify_json_file('output.json', dic_json, time=True)
                break

def annotate_pdb_model_starter(args):
    start_time = time.time()

    dic_json = {'mode':'annotation','parameterization':{'reference':args.pdbmodel, 'strictness_ah':args.strictness_ah,
                'strictness_bs':args.strictness_bs, 'peptide_length':args.peptide_length,
                'width_pic':args.width_pic, 'height_pic':args.height_pic, 'write_graphml':args.write_graphml,
                'write_pdb':True}}
    modify_json_file('output.json', dic_json)


    Bioinformatics.rename_hetatm_and_icode(args.pdbmodel)
    annotate_pdb_model(reference=args.pdbmodel, strictness_ah=args.strictness_ah, strictness_bs=args.strictness_bs,
                       peptide_length=args.peptide_length, width_pic=args.width_pic, height_pic=args.height_pic,
                       write_graphml=args.write_graphml, write_pdb=True, hhr_file = args.hhr_file, verbose=args.verbose,
                       compactness=args.compactness)
    while 1:
        if len(multiprocessing.active_children()) == 0:
            dic_json = {"annotation": "{:.2f}s".format(time.time() - start_time)}
            modify_json_file('output.json', dic_json, time=True)
            break

def decompose_by_community_clustering_starter(args):
    start_time = time.time()

    dic_json = {'mode':'decomposition','parameterization':{'reference':args.pdbmodel, 'strictness_ah':args.strictness_ah,
                                          'strictness_bs':args.strictness_bs, 'peptide_length':args.peptide_length,
                                          'homogeneity':args.homogeneity, 'pack_beta_sheet':args.pack_beta_sheet,
                                          'algorithm':args.algorithm, 'max_ah_dist':args.max_ah_dist,
                                          'min_ah_dist':args.min_ah_dist, 'max_bs_dist':args.max_bs_dist,
                                          'min_bs_dist':args.min_bs_dist, 'write_graphml':args.write_graphml,
                                          'write_pdb':True}}
    modify_json_file('output.json', dic_json)

    Bioinformatics.rename_hetatm_and_icode(args.pdbmodel)
    decompose_by_community_clustering(reference=args.pdbmodel, strictness_ah=args.strictness_ah,
                                      strictness_bs=args.strictness_bs, peptide_length=args.peptide_length,
                                      homogeneity=args.homogeneity,pack_beta_sheet=args.pack_beta_sheet,
                                      max_ah_dist=args.max_ah_dist, min_ah_dist=args.min_ah_dist,
                                      max_bs_dist=args.max_bs_dist, min_bs_dist=args.min_bs_dist,
                                      write_graphml=args.write_graphml, write_pdb=True,  algorithm=args.algorithm,
                                      hhr_file=args.hhr_file, write_graph_img=args.write_graph_img)

    while 1:
        if len(multiprocessing.active_children()) == 0:
            dic_json = {"decomposition": "{:.2f}s".format(time.time() - start_time)}
            modify_json_file('output.json', dic_json, time=True)
            break

def find_local_folds_in_the_graph_starter(args):
    Bioinformatics.rename_hetatm_and_icode(args.pdbmodel)
    find_local_folds_in_the_graph(reference=args.pdbmodel, strictness_ah=args.strictness_ah,
                                  strictness_bs=args.strictness_bs, peptide_length=args.peptide_length, write_pdb=True,
                                  domains=args.domains, compactness=args.compactness, write_graphml=args.write_graphml,
                                  write_graph_img=args.write_graph_img)

def generate_library_starter(args):
    start_time = time.time()

    if args.pdbmodel is not None and os.path.exists(args.pdbmodel):
        Bioinformatics.rename_hetatm_and_icode(args.pdbmodel)
    else:
        args.pdbmodel = None

    if args.ncssearch:
        args.multimer = False

    if not args.score_intra_fragment:
        try:
            args.score_intra_fragment = dict_similarity[args.similarity_intra_chain]
        except:
            args.score_intra_fragment = dict_similarity[3]
    if not args.score_inter_fragments:
        try:
            args.score_inter_fragments = dict_similarity[args.similarity_inter_chains]
        except:
            args.score_inter_fragments = dict_similarity[3]

    dic_json = {'mode': 'library generation', 'submode': 'superposition based on Characteristic Vectors',
                'parameterization': {'reference': args.pdbmodel, 'directory_database':args.directory_database,
                'cath_id':args.cath_id, 'target_sequence':args.target_sequence,'strictness_ah': args.strictness_ah,
                'strictness_bs': args.strictness_bs, 'peptide_length': args.peptide_length,
                'score_intra_fragment':args.score_intra_fragment, 'score_inter_fragments':args.score_inter_fragments,
                'superposition maximum rmsd': args.rmsd_max, 'superposition minimum rmsd':args.rmsd_min,
                'criterium_selection_core':args.criterium_selection_core,
                'clustering_mode':args.clustering_mode, 'number_of_ranges':args.number_of_ranges,
                'number_of_clusters':args.number_of_clusters, 'clustering maximum rmsd':args.rmsd_clustering,
                'local_grid': args.local_grid,'remote_grid': args.remote_grid,'supercomputer': args.supercomputer,
                'use_model_as_it_is':args.use_model_as_it_is}}

    modify_json_file('output.json', dic_json)

    generate_library(local_grid = args.local_grid, remote_grid = args.remote_grid, supercomputer = args.supercomputer, force_core = args.force_core, directory_database = args.directory_database,
    c_angle = args.c_angle, c_dist = args.c_dist, c_angle_dist = args.c_angle_dist, c_cvl_diff = args.c_cvl_diff, score_intra_fragment = args.score_intra_fragment, j_angle = args.j_angle, j_dist = args.j_dist, j_angle_dist = args.j_angle_dist,
    j_cvl_diff = args.j_cvl_diff, score_inter_fragments = args.score_inter_fragments, rmsd_clustering = args.rmsd_clustering, exclude_residues_superpose = args.exclude_residues_superpose, work_directory = args.work_directory,
    targz = args.targz, pdbmodel = args.pdbmodel, remove_coil = args.remove_coil, weight = "distance_avg",
    strictness_ah = args.strictness_ah, strictness_bs = args.strictness_bs, peptide_length = args.peptide_length, enhance_fold = args.enhance_fold, superpose = True, process_join = False,
    nilges = args.nilges, sequence=args.sequence, ncssearch=args.ncssearch, multimer=args.multimer,
    rmsd_min=args.rmsd_min, rmsd_max=args.rmsd_max, ssbridge=args.ssbridge, connectivity=args.connectivity, representative=args.representative, sidechains=args.sidechains,
    criterium_selection_core=args.criterium_selection_core, test=args.test, cath_id=args.cath_id, target_sequence=args.target_sequence,
    clustering_mode=args.clustering_mode, number_of_ranges=args.number_of_ranges, number_of_clusters=args.number_of_clusters, exclude_sequence=args.exclude_sequence,use_model_as_it_is=args.use_model_as_it_is, superposition_mode=False,
    interface=args.interface, verbose=args.verbose)

    while 1:
        if len(multiprocessing.active_children()) == 0:
            dic_json = {"library generation": "{:.2f}s".format(time.time() - start_time)}
            print(dic_json)
            modify_json_file('output.json', dic_json, time=True)
            print()
            if os.path.exists(args.work_directory + "/library/list_rmsd.txt"):
                f = open(args.work_directory + "/library/list_rmsd.txt", "r")
                lines = f.readlines()
                f.close()
            elif os.path.exists(args.work_directory + "/library_cluster/list_rmsd.txt"):
                f = open(args.work_directory + "/library_cluster/list_rmsd.txt", "r")
                lines = f.readlines()
                f.close()
            else:
                print("The file list_rmsd.txt cannot be found, none library could be created")
                raise Exception("The file list_rmsd.txt cannot be found, none library could be created")

            for line in lines:
                line_list = line.split()
                dic_json = {line_list[0]: {'rmsd and core': line_list[1] + "  " + line_list[2]}}
                modify_json_file('output.json', dic_json, lib_generation_superposition=True)

            with open('output.json', 'r') as f:
                dict1 = json.load(f)
            # number_extracted = 0
            # for pdb in dict1['library generation']['extraction']:
            #     number_extracted += dict1['library generation']['extraction'][pdb]
            # mean_extracted_global = round(number_extracted/dict1['library generation']['number of targets'],2)
            # mean_extracted_local = round(number_extracted/len(dict1['library generation']['extraction']),2)
            sum_rmsd = 0.0
            for pdb in dict1['library generation']['superposition']:
                #print(float(dict1['library generation']['superposition'][pdb]['rmsd and core'].split()[0]), type (float(dict1['library generation']['superposition'][pdb]['rmsd and core'].split()[0])))
                sum_rmsd += float(dict1['library generation']['superposition'][pdb]['rmsd and core'].split()[0])
            mean_superposed = round(sum_rmsd/len(dict1['library generation']['superposition']), 2)
            dic_json = { "number of superposed": len(dict1['library generation']['superposition']), "mean of rmsd per superposed": mean_superposed}
            # dic_json = {"number of extracted": number_extracted, "mean of extracted per target":mean_extracted_global,  "mean of extracted per successful extraction":mean_extracted_local,
            #             "number of superposed": len(dict1['library generation']['superposition']), "mean of rmsd per superposed": mean_superposed}
            modify_json_file('output.json', dic_json, lib_generation=True)
            break

def perform_cluster_starter(args):
    """

    :param args:
    """
    cluster_library(directory_database=args.directory_database, rmsd_clustering=args.rmsd_clustering, clustering_mode=args.clustering_mode, number_of_ranges=args.number_of_ranges,
                    number_of_clusters=args.number_of_clusters, exclude_residues_superpose=args.exclude_residues_superpose, ssbridge=args.ssbridge, nilges=args.nilges)

def analyse_interfaces_starter(args):
    start_time = time.time()

    Bioinformatics.rename_hetatm_and_icode(args.pdbmodel)
    analyse_interfaces(reference=args.pdbmodel, strictness_ah=args.strictness_ah, strictness_bs=args.strictness_bs,
                       peptide_length=args.peptide_length, max_ah_dist=args.max_ah_dist, min_ah_dist=args.min_ah_dist,
                       max_bs_dist=args.max_bs_dist, min_bs_dist=args.min_bs_dist, write_graphml=args.write_graphml,
                       write_pdb=True, min_distance=args.min_distance, max_distance=args.max_distance,
                       edges_distance=args.edges_distance, verbose=args.verbose)

def frobenius_starter(args):
    references = []
    if args.reference:
        if os.path.exists(args.reference) and os.path.isfile(args.reference):
            references.append(os.path.abspath(args.reference))
        else:
            print("ERROR: argument reference should be an existing path to a pdb file")
            raise Exception("argument reference should be an existing path to a pdb file")

    elif args.references:
        if os.path.exists(args.references):
            for root, subFolders, files in os.walk(args.references):
                for fileu in files:
                    pdbf = os.path.join(root, fileu)
                    if pdbf.endswith(".pdb"):
                        references.append(os.path.abspath(pdbf))
        else:
            print("ERROR: argument references should be an existing path to a directory containing pdb files")
            raise Exception("argument references should be an existing path to a directory containing pdb files")

    else:
        print("ERROR: either argument reference or argument references should be provided")
        raise Exception("either arguments reference or argument references should be provided")

    if not os.path.isdir(args.targets):
        print("ERROR: argument targets should be an existing path")
        raise Exception("argument targets should be an existing path")
    else:
        targets = [os.path.abspath(os.path.join(args.targets, file)) for file in os.listdir(args.targets) if file.endswith('.pdb')]

    if not os.listdir(args.targets):
        print("ERROR: the directory targets is empty")
        raise Exception("the directory targets is empty")

    if args.distance_min and args.distance_max and args.distance_min>=args.distance_max:
        print("ERROR: argument distance_min should be greater than argument distance_max")
        raise Exception("argument distance_min should be greater than argument distance_max")


    frobenius(references=references, targets=targets, norm_dist=args.distance,
              min_dist=args.distance_min, max_dist=args.distance_max, step_dist=args.distance_step,
              write_plot=args.write_plot, write_matrix=args.write_matrix)

def remove_unconnected_fragments_starter(args):
    remove_unconnected_fragments(reference=args.pdbmodel, strictness_ah=args.strictness_ah,
                                 strictness_bs=args.strictness_bs, peptide_length=args.peptide_length)

#######################################################################################################
#                                             FUNCTIONS                                               #
#######################################################################################################
def merge_dicts(dict1, dict2):
    dict3 = dict1.copy()
    dict3.update(dict2)
    return dict3

def modify_json_file (json_file, dict2, annotation=False, decomposition=False, superposition=False, lib_generation=False,lib_generation_extraction=False, lib_generation_superposition=False, clustering=False, time=False, plot=False, cvs=False):
    if os.path.exists(json_file):
        with open (json_file, 'r') as f:
            dict1 = json.load(f)
    else:
        print("Creating json file")
        with open(json_file, 'w') as f:
            dict1 = {}

    if annotation:
        if 'decomposition' in dict1:
            del dict1['decomposition']
        if 'annotation' in dict1:
            dict1['annotation'] = merge_dicts(dict1['annotation'], dict2) #For python 2 and 3
            #dict1['annotation'] = {**dict1['annotation'], **dict2} #For python 3
        else:
            dict1['annotation'] = dict2
    elif decomposition:
        if 'decomposition' in dict1:
            dict1['decomposition'] = merge_dicts(dict1['decomposition'], dict2)  # For python 2 and 3
            #dict1['decomposition'] = {**dict1['decomposition'], **dict2}
        else:
            dict1['decomposition'] = dict2
    elif lib_generation:
        if not 'library generation' in dict1:
            dict1['library generation'] = {}
            # dict1['library generation']['number of extracted'] = {}
            # dict1['library generation']['mean of extracted per target'] = {}
            # dict1['library generation']['mean of extracted per successful extraction'] = {}
            dict1['library generation']['number of superposed'] = {}
            dict1['library generation']['mean of rmsd per superposed'] = {}
        dict1['library generation'] = merge_dicts(dict1['library generation'],dict2)

    elif lib_generation_extraction:
        if not 'extraction' in dict1['library generation']:
            dict1['library generation']['extraction'] = {}
        dict1['library generation']['extraction'] = merge_dicts(dict1['library generation']['extraction'],dict2)  # For python 2 and 3
    elif lib_generation_superposition:
        if not 'superposition' in dict1['library generation']:
            dict1['library generation']['superposition'] = {}
        dict1['library generation']['superposition'] = merge_dicts(dict1['library generation']['superposition'], dict2)  # For python 2 and 3
        #dict1['library generation']['superposition'] = {**dict1['library generation']['superposition'], **dict2}
    elif clustering:
        if not 'clustering' in dict1['library generation']:
            dict1['library generation']['clustering'] = dict2
        else:
            dict1['library generation']['clustering'] = merge_dicts(dict1['library generation']['clustering'], dict2)  # For python 2 and 3
            #dict1['library generation']['clustering'] = {**dict1['library generation']['clustering'], **dict2}
    elif superposition:
        if 'library generation' in dict1:
            del dict1['library generation']
        if not 'superposition' in dict1:
            dict1['superposition'] = {}
            dict1['superposition'] = dict2
        else:
            dict1['superposition'] = merge_dicts(dict1['superposition'], dict2)  # For python 2 and 3
            #dict1['superposition'] = {**dict1['superposition'], **dict2}
    elif plot:
        if not 'plots' in dict1:
            dict1['plots']={}
            dict1['plots']=dict2
        else:
            dict1['plots'] = merge_dicts(dict1['plots'], dict2)  # For python 2 and 3
            #dict1['plots'] = {**dict1['plots'], **dict2}
    elif cvs:
        if not 'cvs' in dict1:
            dict1['cvs']={}
            dict1['cvs']=dict2
        else:
            dict1['cvs'] = merge_dicts(dict1['cvs'], dict2)  # For python 2 and 3
            #dict1['cvs'] = {**dict1['cvs'], **dict2}
    elif time:
        if not 'time' in dict1:
            dict1['time'] = {}
            dict1['time'] = dict2
        else:
            dict1['time'] = merge_dicts(dict1['time'], dict2)  # For python 2 and 3
            #dict1['time'] = {**dict1['time'], **dict2}
    else:
        dict1 = merge_dicts(dict1, dict2)  # For python 2 and 3
        #dict1 = {**dict1, **dict2}
        # dict1.update(dict2)
        # dict1 = dict(list(dict2.items()) + list(dict1.items()))
    with open(json_file, 'w') as f:
        json.dump(dict1, f)

def galign(s1, s2, gap_open=-10):
    matrix = substitution_matrices.load(name="BLOSUM62")
    gap_open = gap_open  # NOTE #-10
    gap_extend = -1  # NOTE #-0.5
    s1 = s1.replace("-", "X")
    s2 = s2.replace("-", "X")
    s1 = s1.replace("U", "X")
    s2 = s2.replace("U", "X")
    s1 = s1.replace("O", "K")
    s2 = s2.replace("O", "K")
    s1 = s1.replace("B", "X")
    s2 = s2.replace("B", "X")
    s1 = s1.replace("Z", "X")
    s2 = s2.replace("Z", "X")
    s1 = s1.replace("J", "X")
    s2 = s2.replace("J", "X")

    # if "X" in s1 or "X" in s2 or "X" in s3:
    #    return 0

    q1 = pairwise2.align.globalds(s1, s2, matrix, gap_open, gap_extend)

    return q1

def is_compatible(instruction, step, verbose=False, test_equality=False):
    """
    Test if the instruction set and step set correlates

    :param instruction: [continuity:int, (cvl1:float,cvl2:float), angle_between:float, distance:float, angle_distance:float, reslist:list, (sstype1:str,sstype2:str)]
    :param step:        [continuity:int, (cvl1:float,cvl2:float), angle_between:float, distance:float, angle_distance:float, reslist:list, (sstype1:str,sstype2:str)]
    :param thresholds:
    :param verbose:
    :param test_equality:
    :return:
    """

    ins_cont = instruction[0]
    step_cont = step[0]
    uno = numpy.array([instruction[1][0], instruction[1][1]] + instruction[2:5])
    due = numpy.array([step[1][0], step[1][1]] + step[2:5])

    if ins_cont == 1 and step_cont != 1:
        if verbose:
            print("Not Sequential when required", ins_cont, step_cont)
        return False

    corr = scipy.spatial.distance.correlation(uno, due)

    if test_equality:
        thresh = THRESH_CORR[SCALING]["remove_redundance_chain"]
    else:
        thresh = THRESH_CORR[SCALING]["compare_instruction"]

    if verbose:
        print("=====================COMPARISON VECTORS=======================")
        print("uno:", uno)
        print("due:", due)
        print("correlation:", corr, "threshold:", thresh)
        print("==============================================================")

    return corr <= thresh

def compute_instruction(in_a, in_b, angle_type="degree", unique_fragment_cv=False):
    """
    Compute the tertiary relationships between two CVs. Angle, Distance and angle distance

    :param in_a: first CV
    :type in_a: list
    :param in_b: second CV
    :type in_b: list
    :param angle_type: return "degree", "dot", "radians"
    :type angle_type: str
    :return: angle, distance and angle distance
    :rtype: list
    """

    N = numpy.array([[0.0, 0.0, 1.0]])

    _1 = in_a[2] - in_a[3]
    _2 = in_b[2] - in_b[3]

    # Distance Magnitude and distance direction
    D1 = in_b[2] - in_a[2]
    d, D = get_atoms_distance(in_a[2], in_b[2])

    if angle_type.lower() == "degree_by_dot":
        # TETA_1 = angle_between_by_acos(_1, _2)
        TETA_1 = angle_between(_1, _2, N, signed=False)
        TETA_2 = angle_between(_1, D, N, signed=False)
        TETA_2a = angle_between(_2, D, N, signed=False)
        # TETA_2d = angle_between(_1, D1, N, signed=False)
        # TETA_2e = angle_between(_2, D1, N, signed=False)
    elif angle_type.lower() in ["degree", "radians"]:
        TETA_1 = angle_between(_1, _2, N, signed=False)
        if abs(in_a[0] - in_b[0]) == 1:
            TETA_b = angle_between_by_acos(_1, _2)

        TETA_2 = angle_between(_1, D, N, signed=False)
        TETA_2a = angle_between(_2, D, N, signed=False)
        # TETA_2d = angle_between(_1, D1, N, signed=False)
        # TETA_2e = angle_between(_2, D1, N, signed=False)
        # print(TETA_1,TETA_2,TETA_2a,TETA_2d,TETA_2e)
    elif angle_type.lower() == "dot":
        TETA_1 = get_dot_between(_1, _2)
        TETA_2 = TETA_2a = get_dot_between(_1, D)

    TETA_2 = min(TETA_2, TETA_2a)

    if angle_type.lower().startswith("degree"):
        TETA_2 *= 57.2957795
        TETA_1 *= 57.2957795

    if d == 0:
        TETA_2 = 0.0

    valued = []

    if not unique_fragment_cv and in_a[0] == in_b[0]:
        valued = [1.0, (in_a[1], in_b[1]), 0.0, 0.0, 0.0, in_a[4]]
    elif unique_fragment_cv:
        valued = [abs(in_a[0] - in_b[0]), (in_a[1], in_b[1]), TETA_1, d, TETA_2, D, in_a[4] + in_b[4]]
    else:
        valued = [abs(in_a[0] - in_b[0]), (in_a[1], in_b[1]), TETA_1, d, TETA_2, in_a[4] + in_b[4]]

    a = in_a[1]
    b = in_b[1]
    if 2.2 - 0.18 <= a <= 2.2 + 0.18:
        a = "ah"
    elif 1.39 - 0.24 <= a <= 1.39 + 0.24:
        a = "bs"
    else:
        a = "nn"
    if 2.2 - 0.18 <= b <= 2.2 + 0.18:
        b = "ah"
    elif 1.39 - 0.24 <= b <= 1.39 + 0.24:
        b = "bs"
    else:
        b = "nn"

    # valued.append(None)
    valued.append((a, b))
    # if len(in_a) == 5:
    #     in_a.append(None)
    #     in_a.append(a)
    # if len(in_b) == 5:
    #     in_b.append(None)
    #     in_b.append(b)

    return valued


@SystemUtility.timing
def write_image_from_graph(graph, outputname, print_labels=False, x=800, y=800, set_label=None):
    global color_dict

    #NOTE AMB: This code is to customise the colors of the graph
    # list_colors = []
    # dict = {0:'rgb(118,78,128)', 1:'rgb(166,166,166)', 2:'rgb(84, 181, 181)', 3:'rgb(234,86,90)', 4:'rgb(33,140,198)',
    #         5:'rgb(128,128,128)', 6:'rgb(136, 67, 112)', 7:'rgb(50, 86, 139)', 8:'rgb(255,209,35)', 9:'rgb(82, 140, 158)'}
    #
    # for i,clust in enumerate(graph):
    #     for frag in clust:
    #         list_colors.append([frag, dict[i]])
    #
    # new_list_colors =[]
    # list_colors_sorted = sorted(list_colors, key=lambda color: color[0])
    # for num in list_colors_sorted:
    #     new_list_colors.append(num[1])

    if print_labels:
        igraph.plot(graph, outputname, vertex_label=[
            str(vertex["reslist"][0][2]) + "_" + str(vertex["reslist"][0][3][1]) + "-" + str(
                vertex["reslist"][-1][3][1]) for vertex in graph.vs],
                    vertex_color=[color_dict[vertex["reslist"][0][2]] for vertex in graph.vs], vertex_size=60,
                    bbox=(0, 0, x, y))
    elif set_label is not None:
        igraph.plot(graph, outputname, vertex_label=graph.vs[set_label], vertex_size=160, bbox=(0, 0, x, y))
    else:
        igraph.plot(graph, outputname)
        # igraph.plot(graph, outputname, vertex_color = new_list_colors, background=None, vertex_size = 25) #To customise graph

# @SystemUtility.timing
# @SystemUtility.profileit
def format_and_remove_redundance(cvs_global, sep_chains, only_reformat=False):
    """
    Format CVs and remove redundance for chains with identical values of CVs

    :param cvs_global: list of CVs
    :type cvs_global: list of lists
    :param sep_chains: list of delimiters ids for different chains
    :type sep_chains: list of tuples
    :param only_reformat: True if only reformat must be done but not removing, False to allow removing
    :type only_reformat: bool
    :return:
    :rtype:
    :raise: ValueError if sep_chains is empty
    """

    if len(sep_chains) == 0:
        raise ValueError('sep_chains cannot be empty')

    if len(sep_chains) == 1:
        return [cvs_global]

    done = []
    equals = {}
    for i in range(len(sep_chains)):
        a = sep_chains[i]
        if i not in equals:
            equals[i] = []
            done.append(i)
        for j in range(i + 1, len(sep_chains)):
            if j in done:
                continue

            b = sep_chains[j]
            equal = True
            if only_reformat:
                equal = False
            else:
                if a[1] - a[0] == b[1] - b[0]:
                    for p in range(a[1] - a[0] - 1):
                        in_i_a = cvs_global[a[0] + p]
                        in_j_a = cvs_global[a[0] + p + 1]
                        in_i_b = cvs_global[b[0] + p]
                        in_j_b = cvs_global[b[0] + p + 1]
                        dip_a = compute_instruction(in_i_a, in_j_a)
                        dip_b = compute_instruction(in_i_b, in_j_b)

                        if not is_compatible(dip_a, dip_b, verbose=False, test_equality=True):
                            equal = False
                            break
                else:
                    equal = False

            if equal:
                if i in equals:
                    equals[i].append(j)
                else:
                    equals[i] = [j]
                done.append(j)

    cvs_list = []
    for key in equals:
        cvs_list.append(cvs_global[sep_chains[key][0]:sep_chains[key][1]])

    return cvs_list


def validate_residue(residue, ca_list, o_list, all_atom_ca, number_of_residues, ignore):
    """
    Validates a residue and update the passed lists. It allows only valid residues according Bioinformatics.AAList with atoms with at least an occupancy of 10%

    :param residue: Residue object
    :type residue: Bio.PDB.Residue
    :param ca_list: list of coords and properties for CA atoms
    :type ca_list: list [0] pos_X [1] pos_Y [2] pos_Z [3] residue full id [4] residue name code in 3 letters
    :param o_list: list of coords and properties for O atoms
    :type o_list: list [0] pos_X [1] pos_Y [2] pos_Z [3] residue full id [4] residue name code in 3 letters
    :param all_atom_ca: array of CA atoms
    :type all_atom_ca: list of Bio.PDB.Atom
    :param number_of_residues: number of residues read
    :type number_of_residues: int
    :param ignore: set of residue full ids to be ignored
    :type ignore: set
    :return coord_ca, o_list, all_atom_ca, number_of_residues:
    :rtype (list, list, list, int):
    """

    if (residue.get_resname().upper() in Bioinformatics.AAList) and (residue.has_id("CA")) and (
    residue.has_id("O")) and (residue.has_id("C")) and (residue.has_id("N")):
        if residue.get_full_id() in ignore:
            return ca_list, o_list, all_atom_ca, number_of_residues

        ca = residue["CA"]
        o = residue["O"]
        if any(map(lambda x: x.get_occupancy() < 0.1, [ca, o, residue["C"], residue["N"]])):
            return ca_list, o_list, all_atom_ca, number_of_residues

        number_of_residues += 1

        # for each atom 4 values are saved:
        # [0] pos_X [1] pos_Y [2] pos_Z [3] residue full id [4] residue name code in 3 letters
        co_ca = ca.get_coord()
        co_o = o.get_coord()
        ca_list.append(
            [float(co_ca[0]), float(co_ca[1]), float(co_ca[2]), residue.get_full_id(), residue.get_resname()])
        all_atom_ca.append(ca)
        o_list.append([float(co_o[0]), float(co_o[1]), float(co_o[2]), residue.get_full_id(), residue.get_resname()])


    return ca_list, o_list, all_atom_ca, number_of_residues


def get_cvs(structure, use_list=[], ignore_list=[], length_fragment=3, one_model_per_nmr=False, process_only_chains=[]):
    """
    Parse a structure and generate for every consecutive 'length_fragment' residues a CV with a window f 1 overlapping residue.

    :param structure: The structure object of the structure
    :type structure: Bio.PDB.Structure
    :param use_list:
    :type use_list: list
    :param ignore_list:
    :type ignore_list: list
    :param length_fragment: Define the peptide length for computing a CV
    :type length_fragment: int
    :param one_model_per_nmr:
    :type one_model_per_nmr: bool
    :return (cvs, separating_chains): List of CVs and the id separation for each chain
    :rtype: (list,list)
    """

    allAtomCA = []
    ca_list = []
    o_list = []
    numberOfResidues = 0
    separating_chains = []
    start_separa = 0
    end_separa = 0
    cvs = []

    if len(use_list) > 0:
        for fra in use_list:
            lir = [residue for residue in Bio.PDB.Selection.unfold_entities(structure[fra["model"]][fra["chain"]], "R")
                   if residue.get_id() in fra["resIdList"]]
            for residue in lir:
                if process_only_chains is not None and len(process_only_chains) > 0 and residue.get_full_id()[2] not in process_only_chains: continue
                ca_list, o_list, allAtomCA, numberOfResidues = validate_residue(residue, ca_list, o_list, allAtomCA,
                                                                                numberOfResidues, ignore_list)
    elif one_model_per_nmr:
        #print('SHERLOCK structure.get_list()[0]',structure.get_list()[0])
        #print('SHERLOCK type(structure.get_list()[0])', type(structure.get_list()[0]))
        lir = Bio.PDB.Selection.unfold_entities(structure.get_list()[0], "R")
        for residue in lir:
            if process_only_chains is not None and len(process_only_chains) > 0 and residue.get_full_id()[2] not in process_only_chains: continue

            ca_list, o_list, allAtomCA, numberOfResidues = validate_residue(residue, ca_list, o_list, allAtomCA,
                                                                            numberOfResidues, ignore_list)
    else:
        lir = Bio.PDB.Selection.unfold_entities(structure, "R")
        for residue in lir:
            if process_only_chains is not None and len(process_only_chains) > 0 and residue.get_full_id()[2] not in process_only_chains: continue

            ca_list, o_list, allAtomCA, numberOfResidues = validate_residue(residue, ca_list, o_list, allAtomCA,
                                                                            numberOfResidues, ignore_list)
    if numberOfResidues < 3:
        print("ERROR: There are not enough residues with occupancy higher than the 10% for computing Characteristic Vectors")
        raise Exception("ERROR: There are not enough residues with occupancy higher than the 10% for computing Characteristic Vectors")
    # This is the number of possible cvs n_cvs = n_aa - (n_peptides - 1) = n_aa - n_peptides + 1
    numberOfSegments = numberOfResidues - length_fragment + 1

    coordCA = numpy.array([c[:3] for c in ca_list])
    coordO = numpy.array([c[:3] for c in o_list])

    if numberOfSegments <= 0:
        # print ("\t\t\tNo enough residues available to create a fragment")
        return cvs, separating_chains

    vectorsCA = numpy.empty((numberOfSegments, 3))
    vectorsCA.fill(0.0)
    vectorsO = numpy.empty((numberOfSegments, 3))
    vectorsO.fill(0.0)
    vectorsH = numpy.empty((numberOfSegments, 1))
    vectorsH.fill(0.0)

    # a  b  c  d  e  f
    # |--0--|
    #    |--1--|
    #       |--2--|
    #          |--3--|
    # ==================
    # 4 vectors (from 0 to 3)
    #
    # |--0--| = a/3+b/3+c/3
    # |--1--| = |--0--|+(d-a)/3 = a/3 + b/3 + c/3 + d/3 -a/3 = b/3 + c/3 + d/3
    # |--2--| = |--1--|+(e-b)/3 = b/3 + c/3 + d/3 + e/3 -b/3 = c/3 + d/3 + e/3

    vectorsCA[0] = vectorsCA[0] + (coordCA[:length_fragment, :] / float(length_fragment)).sum(axis=0)
    vectorsO[0] = vectorsO[0] + (coordO[:length_fragment, :] / float(length_fragment)).sum(axis=0)
    for i in range(1, len(vectorsCA)):
        vectorsCA[i] = vectorsCA[i - 1] + (coordCA[i + length_fragment - 1] - coordCA[i - 1]) / float(length_fragment)
        vectorsO[i] = vectorsO[i - 1] + (coordO[i + length_fragment - 1] - coordO[i - 1]) / float(length_fragment)

    H = vectorsCA - vectorsO
    vectorsH = numpy.sqrt((H ** 2).sum(axis=1))

    last_chain = None

    for i in range(len(vectorsCA)):
        # if same occupancy take the lowest bfactor
        prevRes = (" ", None, " ")
        ncontigRes = 0
        resids = []
        prev_model = None
        prev_chain = None

        for yui in range(i, i + length_fragment):  # quindi arrivo a i+lengthFragment-1
            resan = (ca_list[yui])[3]
            resids.append(list(resan) + [ca_list[yui][4]])
            resa = resan[3]

            if prevRes == (" ", None, " "):
                ncontigRes += 1
            elif prev_chain is None or prev_chain == resan[2]:
                resaN = Bioinformatics.get_residue(structure, resan[1], resan[2], resan[3])
                prevResC = Bioinformatics.get_residue(structure, prev_model, prev_chain, prevRes)
                if Bioinformatics.check_continuity(prevResC, resaN):
                    ncontigRes += 1
            prevRes = resa
            prev_model = resan[1]
            prev_chain = resan[2]

        if ncontigRes != length_fragment:
            vectorsH[i] = 100  # this value identify a not reliable measure for cv
        else:
            cvs.append([i, vectorsH[i], vectorsCA[i], vectorsO[i], resids])

        if prev_chain == last_chain:
            end_separa = len(cvs)
        else:
            last_chain = prev_chain
            if start_separa < end_separa:
                separating_chains.append((start_separa, end_separa))
            start_separa = end_separa

    if start_separa < end_separa:
        separating_chains.append((start_separa, end_separa))

    if len(separating_chains) > 0 and separating_chains[0] == (0, 0):
        separating_chains = separating_chains[1:]  # this is done to eliminate the first (0,0) that is unuseful

    if len(cvs) == 0:
        print("ERROR: Characteristics Vectors cannot be computed.")
        raise Exception("Characteristics Vectors cannot be computed")

    return cvs, separating_chains


def get_ss_from_cvl(cvl1):
    """
    Return ah,bs or coil for the cvl in input

    :param cvl1:
    :type cvl1:
    :return:
    :rtype:
    """

    delta_cvla = 0.2  # 0.2
    delta_cvlb = 0.05  # 0.12
    exty1 = "bs" if (1.4 - delta_cvlb <= float(cvl1) <= 1.4 + delta_cvlb) else "coil"
    if exty1 == "coil":
        exty1 = "ah" if (2.2 - delta_cvla <= float(cvl1) <= 2.2 + delta_cvla) else "coil"
    return exty1


def get_unique_cv_among_residues(strucc, resilist):
    """

    :param strucc:
    :type strucc:
    :param resilist:
    :type resilist:
    :return:
    :rtype:
    """

    listCA = numpy.empty((len(resilist), 3))
    listO = numpy.empty((len(resilist), 3))
    for i, x in enumerate(resilist):
        resi = Bioinformatics.get_residue(strucc, x[1], x[2], x[3])
        listCA[i] = resi["CA"].get_coord()
        listO[i] = resi["O"].get_coord()
    if len(listCA) == 0 or len(listO) == 0:
        return None

    CAx = listCA.mean(axis=0)
    Ox = listO.mean(axis=0)
    return [0, 0, CAx, Ox, []]


#def get_all_fragments(graph: igraph.Graph) -> typing.List[igraph.Vertex]:
def get_all_fragments(graph):

    graph.vs["resIdList"] = graph.vs["reslist"]
    graph.vs["pdbid"] = [p[0][0] for p in graph.vs["reslist"]]
    graph.vs["model"] = [p[0][1] for p in graph.vs["reslist"]]
    graph.vs["chain"] = [p[0][2] for p in graph.vs["reslist"]]
    graph.vs["fragLength"] = [len(p) for p in graph.vs["reslist"]]
    graph.vs["index"] = [v.index for v in graph.vs]
    l = sorted([v.attributes() for v in graph.vs], key=lambda x: tuple(x["reslist"][0]))

    return l


def get_connected_fragment_to_edge(fragment, edge):
    ind = fragment.index
    if edge.source == ind:
        return edge.graph.vs[edge.target]
    elif edge.target == ind:
        return edge.graph.vs[edge.source]
    else:
        return None


def get_vseq_neighbours_fragments(graph, fragment, sortmode=None):
    if not sortmode:
        return fragment.neighbors()  # VertexSeq
    else:
        return [get_connected_fragment_to_edge(fragment, edge) for edge in
                sorted(graph.es.select(graph.incident(fragment)), key=lambda x: x[sortmode])]  # VertexSeq


def get_eseq_neighbours_fragments(graph, fragment, sortmode=None):  # sortmode="avg",sortmode="min"
    if not sortmode:
        return [graph.es[graph.get_eid(frag.index, fragment.index)] for frag in fragment.neighbors()]  # EdgeSeq
    else:
        #print(graph.incident(fragment), type(graph.es.select(graph.incident(fragment))))
        return sorted(graph.es.select(graph.incident(fragment)), key=lambda x: x[sortmode])  # EdgeSeq


@SystemUtility.timing
# @SystemUtility.profileit
def get_3d_cvs_matrix(cvs, is_model, maximum_distance=None, maximum_distance_bs=None, just_diagonal_plus_one=False,
                      mixed_chains=False):
    if mixed_chains:
        cvs = [cv for cvst in cvs for cv in cvst]

    n = len(cvs)

    # matrice = ADT.get_matrix(n,n)
    high_d = 0

    ####Sparse matrices we just need a sparse triangular matrix that can grow
    ###matrice = scipy.sparse.lil_matrix((n,n), dtype=list)

    #TODO: Please pay extremely attention if just_diagonal_plus=True this code will not work because of the new annotation system
    #TODO: needs to explore several regions of the matrix to define the tertiary structure porperties.

    if is_model:
        matrix = {}
        matrix["n"] = n
        matrix["fragments_bounds"] = [[0]]

        for i in range(n):
            for j in range(i, n):
                if just_diagonal_plus_one and j != i + 1:
                    continue

                valued = compute_instruction(cvs[i], cvs[j])
                if maximum_distance is not None and valued[3] >= maximum_distance:
                    continue

                if maximum_distance_bs is not None and valued[-1] == ("bs", "bs") and valued[3] >= maximum_distance_bs:
                    continue

                matrix[(i, j)] = valued
                if valued[3] > high_d:
                    high_d = valued[3]

                if j == i + 1:
                    if valued[0] != 1:
                        matrix["fragments_bounds"][-1].append(i) #TODO: Check if I am using the indeces of the cvs_list array or the values because if it is the second this should be cvs[i]
                        matrix["fragments_bounds"].append([j])
    else:
        # NOTE: when is not is_model maximum_distance and maximum_distance_bs are ignored. If we really want that filters then after this creation
        # I should do something to filter them out.
        matrix = {(i, j): compute_instruction(cvs[i], cvs[j]) for i in range(n) for j in range(i, n) if
                  not just_diagonal_plus_one or j == i + 1}
        matrix["n"] = n

    if is_model:
        if matrix[(n - 2, n - 1)][0] == 1:
            matrix["fragments_bounds"][-1].append(n - 1)
        else:
            matrix["fragments_bounds"][-1].append(n - 1)
            # matrice["fragments_bounds"].append([n-1,n-1])
        return matrix, cvs, high_d
    else:
        return matrix, cvs, None

def print_pattern(pattern):
    print("CONTINOUS FRAGMENTS:")
    for f in range(len(pattern["fragments_bounds"])):
        print("=====================Frag. n.: ", f, "========================")
        index = pattern["fragments_bounds"][f]
        for t in range(index[0],index[1]):
            i, j = t, t+1
            u = pattern[(i, j)]

            if i != j:
                ssdes = ""
                if u[6][0] == u[6][1] and u[0] == 1:
                    ssdes = u[6][0] + "\t"
                else:
                    ssdes = u[6][0] + " - " + u[6][1] + "\t"

                print(ssdes, i, j, u[0], u[1], u[2], u[3], u[4], u[5][0][1], u[5][0][2], list(map(lambda x: x[3][1], u[5])))

        print("==============================================================")

    print()

@SystemUtility.timing
# @SystemUtility.profileit
def aleph_secstr(strucc, cvs_list, matrix, min_ah=None, min_bs=None, strictness_ah=0.45, strictness_bs=0.20, interface = False):
    """
    Annotates the secondary structure through CVs and return a graph with no edges
    # FORMULA ALPHA ANGLE
    #n = 0.5
    #mean = 20.0
    #v = 0.9
    #p = 0.5
    #topx = 180.0
    #step = 0.01

    # FORMULA BETA ANGLE
    #n=0.5
    # mean = 54.0
    # v = 0.9
    # p = 0.5
    # topx = 180.0
    # step = 0.01

    # FORMULA ALPHA CVL
    #n=0.5 or 1
    # mean = 2.2
    # v = 0.9
    # v = 0.0 this is to have a plateau with 0 after 2.2
    # p = 0.8
    # p = 1.0 this is to have a plateau with 0 after 2.2
    # topx = 2.4
    # step = 0.01

    # FORMULA BETA CVL
    # n=0.5 or 1
    # mean = 1.4
    # v = 0.9
    # p = 0.8
    # topx = 2.4
    # step = 0.01

    :param strucc:
    :type strucc:
    :param cvs_list:
    :type cvs_list:
    :param min_ah:
    :type min_ah:
    :param min_bs:
    :type min_bs:
    :param strictness_ah:
    :type strictness_ah:
    :param strictness_bs:
    :type strictness_bs:
    :return:
    :rtype:
    """

    global BS_UD_EA
    global BS_UU_EA
    global BS_MAX

    dist_mean_bs = 5.1
    dist_mean_ah = 0.0
    dist_num_ah = 0
    angle_mean_bs = 54
    angle_mean_ah = 20
    cvl_mean_bs = 1.4
    cvl_mean_ah = 2.2
    thresh_scores = 1.5

    a = [[lis[1], get_ss_from_cvl(lis[1]), lis[4], lis[0]] for lis in cvs_list]
    dimap = {value[0]: i for (i, value) in enumerate(cvs_list)}

    def __scoring_tick_fn(u, mean, topx, v=0.9, p=0.5, n=0.5):
        r = numpy.abs((u - mean) / mean) ** n if u <= mean else (((((u - mean) * ((mean * v))) / (
        mean + ((topx - mean) * p) - mean)) / mean)) ** n if u <= mean + ((topx - mean) * p) else (v + (((((u - (
        mean + (topx - mean) * p)) * (mean - (mean * v))) / (topx - (mean + (topx - mean) * p) + (
        mean * v))) / mean))) ** n
        return r

    def __check_by_unified_score_step2(uno, due, dizio3d, take_first=True, validate=["bs", "coil"], min_num_bs=1.0):
        ###value1 = compute_instruction(cvs_list[dimap[uno[3]]], cvs_list[dimap[due[3]]])
        sup = tuple(sorted([dimap[uno[3]], dimap[due[3]]]))
        value1 = matrix[sup]
        if value1[0] != 1:
            return uno, due, True

        beta_score_uno = 100000
        alpha_score_uno = 100000
        beta_score_due = 100000
        alpha_score_due = 100000
        #if uno[3] == 134:
        #    print("BEFORE UNO",uno[1], uno[0], value1[2], uno[2], "---", uno[3])
        if take_first and (dizio3d is None or uno[1] in validate):
            beta_score_uno = __scoring_tick_fn(uno[0], cvl_mean_bs, 2.4, v=0.9, p=0.8, n=1) + __scoring_tick_fn(
                value1[2], angle_mean_bs, 180.0, v=0.9, p=0.5, n=0.5) #n=0.5
            alpha_score_uno = __scoring_tick_fn(uno[0], cvl_mean_ah, 2.4, v=0.0, p=1.0, n=1) + __scoring_tick_fn(
                value1[2], angle_mean_ah, 180.0, v=0.9, p=0.5, n=0.5)
            CA_CA_d = value1[3]

            if dizio3d is not None:
                #if uno[3] == 134:
                #    print("UNO cvl+int_angles bs:",beta_score_uno,"ah:",alpha_score_uno)
                if uno[3] not in dizio3d or len(dizio3d[uno[3]]) == 0:
                    listuno = []
                    beta_score_uno += 5
                else:
                    listuno = sorted(dizio3d[uno[3]], key=lambda x: x[2])
                    # NOTE: Here we do not take the absolute value of the difference because when the min distance is lower than the mean
                    #      in principle is not a bad condition and we want to add just 0 in that case
                    beta_score_uno += __scoring_tick_fn(listuno[0][2], dist_mean_bs, 10.0, v=0.9, p=0.5, n=0.5) if listuno[0][2] >= dist_mean_bs else 0  #n=0.5

                if len(listuno) < min_num_bs:
                    beta_score_uno += 5
                else:
                    # NOTE: Add the external angles
                    listuno = sorted(listuno, key=lambda x: x[2], reverse=True)
                    listunobleah = sorted([p[2] for p in listuno], reverse=True)
                    f = [((-1.0 * (t + 1)) / (e * len(listuno)))*2*(max(BS_UD_EA[int(round(listuno[t][1]))],BS_UU_EA[int(round(listuno[t][1]))])/BS_MAX[numpy.argmax([BS_UD_EA[int(round(listuno[t][1]))],BS_UU_EA[int(round(listuno[t][1]))]])]) for t, e in enumerate(listunobleah)]
                    #print("f is", f)
                    #print("Ho sottratto sum(f)", 1 + sum(f), f)
                    beta_score_uno += sum(f)

                beta_score_uno /= 4
                alpha_score_uno /= 2
                #if uno[3] == 134:
                #    print("UNO cvl+int_angles+dist+numlinks bs:",beta_score_uno,"ah:",alpha_score_uno,"len(listuno)",len(listuno),"f",sum(f))
                ###print("UNO: NUMLINKS BS",len(listuno))

            if len(uno) == 4:
                uno.append(numpy.abs(alpha_score_uno - beta_score_uno))
            else:
                uno[-1] = numpy.abs(alpha_score_uno - beta_score_uno)

            if beta_score_uno < alpha_score_uno and beta_score_uno <= thresh_scores and numpy.abs(alpha_score_uno - beta_score_uno) >= strictness_bs:
                uno[1] = "bs"
            elif dizio3d is None and alpha_score_uno < beta_score_uno and alpha_score_uno <= thresh_scores and numpy.abs(
                            alpha_score_uno - beta_score_uno) >= strictness_ah:
                uno[1] = "ah"
            elif dizio3d is not None and alpha_score_uno < beta_score_uno and alpha_score_uno <= thresh_scores and numpy.abs(
                            alpha_score_uno - beta_score_uno) >= strictness_ah:
                uno[1] = "coil"
            else:
                uno[1] = "coil"

            #if uno[3] == 134:
            #    print("ANNOTATION UNO IS:",uno)

        if (dizio3d is None or due[1] in validate):
            ###print("CVL score:",__scoring_tick_fn(due[0],cvl_mean_bs,2.4,v=0.9,p=0.8,n=1),"con cvl=",due[0])
            ###print("Angle score:",__scoring_tick_fn(value1[2],angle_mean_bs,180.0,v=0.9,p=0.5,n=0.5),"con angle=",value1[2])
            beta_score_due = __scoring_tick_fn(due[0], cvl_mean_bs, 2.4, v=0.9, p=0.8, n=1) + __scoring_tick_fn(
                value1[2], angle_mean_bs, 180.0, v=0.9, p=0.5, n=0.5) #n=0.5
            alpha_score_due = __scoring_tick_fn(due[0], cvl_mean_ah, 2.4, v=0.0, p=1.0, n=1) + __scoring_tick_fn(
                value1[2], angle_mean_ah, 180.0, v=0.9, p=0.5, n=0.5)
            CA_CA_d = value1[3]
            ###print("XCa-XCa distance",CA_CA_d)
            #if due[3] == 134:
            #    print("BEFORE DUE", due[1], due[0], value1[2], due[2], "---", due[3])
            #if due[3] == 134:
            #    print("DUE cvl+int_angles bs:", beta_score_due, "ah:", alpha_score_due)
            #    print("CVL_BS",__scoring_tick_fn(due[0], cvl_mean_bs, 2.4, v=0.9, p=0.8, n=1))
            #    print("CVL_AH",__scoring_tick_fn(due[0], cvl_mean_ah, 2.4, v=0.0, p=1.0, n=1))
            #    print("ANG_BS",__scoring_tick_fn(value1[2], angle_mean_bs, 180.0, v=0.9, p=0.5, n=0.5))
            #    print("ANG_AH",__scoring_tick_fn(value1[2], angle_mean_ah, 180.0, v=0.9, p=0.5, n=0.5))

            if dizio3d is not None:
                ###print(due,"VALIDATE is",validate)
                ###print("DUE BEFORE bs:", beta_score_due, "ah:", alpha_score_due)
                if due[3] not in dizio3d or len(dizio3d[due[3]]) == 0:
                    #print("We are HEre where list is 0")
                    listdue = []
                else:
                    listdue = sorted(dizio3d[due[3]], key=lambda x: x[2])
                    ###print("How is it listdue",[(q[0],q[1]) for q in sorted(dizio3d[due[3]], key=lambda x: x[0])])
                    # NOTE: Here we do not take the absolute value of the difference because when the min distance is lower than the mean
                    #      in principle is not a bad condition and we want to add just 0 to the score
                    beta_score_due += __scoring_tick_fn(listdue[0][2], dist_mean_bs, 10.0, v=0.9, p=0.5, n=0.5) if \
                    listdue[0][2] >= dist_mean_bs else 0 #n=0.5

                if len(listdue) < min_num_bs:
                    #print("MIN is",min_num_bs,"LEN is",len(listdue),listdue,"BETASCORE was",beta_score_due)
                    beta_score_due += 5
                else:
                    #TODO: Add the external angles
                    listdue = sorted(listdue, key= lambda x: x[2], reverse=True)
                    #print(listdue)
                    # for t,e in enumerate(listdue):
                    #     y = (-1.0 * (t + 1)) / (e[2] * len(listdue))
                    #     print("Base y=",y)
                    #     m = max(BS_UD_EA[int(round(listdue[t][1]))],BS_UU_EA[int(round(listdue[t][1]))])
                    #     print("The angle is",int(round(listdue[t][1])),"in BS_UD f is",BS_UD_EA[int(round(listdue[t][1]))],"in BS_UU f is",BS_UU_EA[int(round(listdue[t][1]))],"max is",m)
                    #     f = BS_MAX[numpy.argmax([BS_UD_EA[int(round(listdue[t][1]))],BS_UU_EA[int(round(listdue[t][1]))]])]
                    #     print("MAX of the corresponding distribution ",f)
                    #     print("Value",y*(m/f))
                    f = [((-1.0 * (t + 1)) / (e[2] * len(listdue)))*2*(max(BS_UD_EA[int(round(listdue[t][1]))],BS_UU_EA[int(round(listdue[t][1]))])/BS_MAX[numpy.argmax([BS_UD_EA[int(round(listdue[t][1]))],BS_UU_EA[int(round(listdue[t][1]))]])]) for t, e in enumerate(listdue)]
                    #print("BETA SCORE WAS:",beta_score_due)
                    beta_score_due += sum(f)
                    #print("f is", f)
                    #print("Ho sottratto sum(f)", 1 + sum(f), f)
                    #print("Ho sottratto sum(f)",sum(f), f)
                ###cont_dist = [p[0] for p in sorted(listdue, key=lambda x: x[0])]
                ###num_cont_dist = sum([0 if o == 0 or cont_dist[o - 1] + 1 == u else 1 for o, u in enumerate(cont_dist)])
                ###print("Nuovo parametro continuita distanze",cont_dist,num_cont_dist)

                beta_score_due /= 4
                alpha_score_due /= 2

                #if due[3] == 134:
                #    print("DUE cvl+int_angles+dist+numlinks bs:",beta_score_uno,"ah:",alpha_score_uno,"len(listuno)",len(listdue),"f",sum(f))
                ###print("DUE: NUMLINKS BS",len(listdue))

            if len(due) == 4:
                due.append(numpy.abs(alpha_score_due - beta_score_due))
            else:
                due[-1] = numpy.abs(alpha_score_due - beta_score_due)

            if beta_score_due < alpha_score_due and beta_score_due <= thresh_scores and numpy.abs(
                            alpha_score_due - beta_score_due) >= strictness_bs:
                due[1] = "bs"
            elif dizio3d is None and alpha_score_due < beta_score_due and alpha_score_due <= thresh_scores and numpy.abs(
                            alpha_score_due - beta_score_due) >= strictness_ah:
                due[1] = "ah"
            elif dizio3d is not None and alpha_score_due < beta_score_due and alpha_score_due <= thresh_scores and numpy.abs(
                            alpha_score_due - beta_score_due) >= strictness_ah:
                due[1] = "coil"
            else:
                due[1] = "coil"

        # print("UNO AFTER bs:", beta_score_uno, "ah:", alpha_score_uno)
        # print(uno)
        # print()
        # #
        # print("DUE AFTER bs:", beta_score_due, "ah:", alpha_score_due)
        # print(due)
        # print()
                #if due[3] == 134:
            #    print("ANNOTATION DUE IS:",due)

        return uno, due, False

    def __get_associations(a, stringent=True):
        associations = {}
        for e in a:
            for resi in e[2]:
                if tuple(resi) not in associations:
                    associations[tuple(resi)] = {"bs": 0, "ah": 0, "coil": 0, "COIL": 0}
                associations[tuple(resi)][e[1]] += 1

        lisorted = sorted(associations.keys(), key=lambda x: x[:3]+(x[3][1:],))
        for t, key in enumerate(lisorted):
            result = ""
            # or (associations[key]["ah"] >= 2 and associations[key]["bs"] == 0) \
            # or (t>0 and t<len(associations.keys())-1 and associations[key]["ah"] == 1 and associations[key]["coil"] == 2 and associations[lisorted[t-1]]["result"] == "coil" and associations[lisorted[t+1]]["ah"] >= 2 and associations[lisorted[t+1]]["bs"] == 0):
            max_for_key = sum([associations[key][z] for z in associations[key].keys()])
            if max_for_key > 3: max_for_key = 3

            if associations[key]["ah"] == max_for_key \
                    or (max_for_key < 3 and  associations[key]["ah"] >=1 and associations[key]["bs"] == 0) or (
                    0 < t < len(associations.keys()) - 1 and associations[key]["ah"] >= 2 and
                    associations[key]["bs"] == 0 and associations[lisorted[t - 1]]["result"] == "coil" and
                    associations[lisorted[t + 1]]["ah"] >= 3):
                result = "ah"
            elif stringent and associations[key]["bs"] == max_for_key:
                result = "bs"
            elif not stringent and (associations[key]["bs"] == max_for_key
                                    or (associations[key]["bs"] == 2 and associations[key]["COIL"] == 0)
                                    or (associations[key]["bs"] == 1 and associations[key]["ah"] < 2 and associations[key]["COIL"] == 0)
                                    or (associations[key]["bs"] >= 1 and associations[key]["COIL"] == 1 and t > 0 and
                                        associations[key]["bs"] >= 1 >= associations[key]["ah"] and associations[lisorted[t - 1]]["result"] in ["coil", "COIL"])
                                    or (associations[key]["bs"] >= 1 and associations[key]["COIL"] == 1 and t > 1 and associations[key]["bs"] >= 1 and associations[key]["ah"] <= 1 and associations[lisorted[t - 2]]["result"] in ["coil", "COIL"])
                                    or (associations[key]["bs"] >= 1 and associations[key]["COIL"] == 1 and t < len(associations.keys()) - 1 and associations[key]["bs"] >= 1 and
                                        associations[key]["ah"] <= 1 < associations[lisorted[t + 1]]["coil"] +
                                        associations[lisorted[t + 1]]["COIL"])
                                    or (associations[key]["bs"] >= 1 and associations[key]["COIL"] == 1 and t < len(associations.keys()) - 2 and associations[key]["bs"] >= 1 and
                                        associations[key]["ah"] <= 1 < associations[lisorted[t + 2]]["coil"] +
                                        associations[lisorted[t + 2]]["COIL"])):

                # or (t>0 and t<len(associations.keys())-1 and associations[key]["bs"] == 1 and associations[key]["coil"] == 2 and associations[lisorted[t-1]]["result"] == "coil" and associations[lisorted[t+1]]["bs"] >= 2 and associations[lisorted[t+1]]["ah"] == 0)):
                result = "bs"
            else:
                result = "coil"

            #print("STRINGENT:", stringent, "T all:", associations[key], "T:", t, "RESULT:",result,"RESI",key)  # ,"T-1 all: "+str(associations[lisorted[t - 1]]) if t>0 else "","T+1 all: "+str(associations[lisorted[t+1]]) if t<len(associations.keys()) - 1 else "")
            associations[key]["result"] = result
            #print(key, associations[key], result)
        return associations

    def __generate_graph2(a, min_ah, min_bs, associations):
        if min_ah is not None and min_ah < 3:
            min_ah = 3
        if min_bs is not None and min_bs < 3:
            min_bs = 3

        listaFrags = []
        used = set([])
        for e in a:
            #print(e)
            if len(listaFrags) == 0:
                #print("A")
                listaFrags.append([])
                #print("Was empty adding []",listaFrags)
            if e[1].lower() != "coil" and len(listaFrags[-1]) == 0:
                #print("B")
                #print("e[1]=",e[1],"that is not coil and ",len(listaFrags[-1]),"is empty")
                listaFrags[-1].append(e)
                for r in e[2]:
                    used.add(tuple(r))
            elif any(map(lambda x: x[1].lower() == e[1].lower() and x[2][0][2] == e[2][0][2] and len(
                    (set([tuple(p) for p in x[2]]) & set([tuple(p) for p in e[2]]))) > 0, listaFrags[-1])):
                #print("C")
                #print("e[1]=",e[1],"but could be automatically attachable to ",listaFrags[-1])
                listaFrags[-1].append(e)
                for r in e[2]:
                    used.add(tuple(r))
            elif any(map(lambda x: e[1].lower() != "coil" and x[1] in ["ah", "bs"] and x[1] != e[1] and x[2][0][2] == e[2][0][
                2] and len((set([tuple(p) for p in x[2]]) & set([tuple(p) for p in e[2]]))) > 0, listaFrags[-1])):
                #print("We need to change", e[1], "into coil because previous was a different sstype and it is continous")
                #print("e[1]=",e[1],"that is not coil and we need to open a new list")
                e[1] = "coil"
                listaFrags.append([e])
                for r in e[2]:
                    used.add(tuple(r))
            elif e[1].lower() != "coil":  # and all(map(lambda x: tuple(x) not in used, e[2])):
                #print("D")
                #print("e[1]=",e[1],"that is not coil and we need to open a new list")
                listaFrags.append([e])
                for r in e[2]:
                    used.add(tuple(r))
            elif e[1].lower() == "coil":
                 #print (e[2],"this is coil")
                 listaFrags.append([e])
                 for r in e[2]:
                     used.add(tuple(r))
            else:
                #print ("E")
                #print ("e[1]=", e[1], "creating an empty list")
                listaFrags.append([e])

        # print("===============================LISTAFRAGS===============================")
        # for fra in listaFrags:
        #     print(fra)
        # print("=========================================================================")

        listaFrags = [
            {"sstype": fragl[0][1].lower(), "reslist": fragl[0][2][:] + list(map(lambda x: x[2][-1], fragl[1:])),
             "cvids": list(map(lambda x: x[3], fragl)), "cvls": list(map(lambda x: x[0], fragl)), "sequence": "".join(
                list(map(lambda x: Bioinformatics.AADICMAP[x[4]],
                         fragl[0][2][:] + list(map(lambda x: x[2][-1], fragl[1:])))))} for
            fragl in listaFrags if len(fragl) > 0 and fragl[0][1] in ["ah", "bs", "coil","COIL"]]

        #
        # for frag in listaFrags:
        #         print (frag["sstype"])

        seen = set([])
        for o, frag in enumerate(listaFrags):
            newcorrect = []
            newseq = ""
            resforgraph = []
            seqforgraph = ""
            for resi in frag["reslist"]:
                #print (resi, associations[tuple(resi)]["result"], "==", frag["sstype"])
                if associations[tuple(resi)]["result"] == frag["sstype"] and tuple(resi) not in seen:
                    if o > 0 and len(listaFrags[o - 1]["reslist"]) > 0:
                        # print (listaFrags[o-1])
                        resaN = Bioinformatics.get_residue(strucc, resi[1], resi[2], resi[3])
                        prev = listaFrags[o - 1]["reslist"][-1]
                        prevResC = Bioinformatics.get_residue(strucc, prev[1], prev[2], prev[3])
                        if Bioinformatics.check_continuity(resaN, prevResC) and listaFrags[o - 1]["sstype"] in ["ah", "bs"]:
                            resforgraph.append((0, "-", "-", ('', '-', '')))
                            seqforgraph += "-"
                        else:
                            newcorrect.append(resi)
                            newseq += Bioinformatics.AADICMAP[resi[4]]
                            resforgraph.append(resi)
                            seqforgraph += Bioinformatics.AADICMAP[resi[4]]
                            if tuple(resi) not in seen:
                                seen.add(tuple(resi))
                    else:
                        newcorrect.append(resi)
                        newseq += Bioinformatics.AADICMAP[resi[4]]
                        resforgraph.append(resi)
                        seqforgraph += Bioinformatics.AADICMAP[resi[4]]
                        if tuple(resi) not in seen:
                            seen.add(tuple(resi))
                else:
                    resforgraph.append((0, "-", "-", ('', '-', '')))
                    seqforgraph += "-"
            frag["reslist"] = newcorrect
            frag["sequence"] = newseq
            frag["resforgraph"] = resforgraph
            frag["seqforgraph"] = seqforgraph

            if len(frag["reslist"]) == 0:
                frag["sstype"] = "coil"

        if min_ah is not None and min_bs is not None:
            listaFrags = [fr for i, fr in enumerate(listaFrags) if
                          (fr["sstype"] == "coil" and len(fr["reslist"]) > 0) or
                          (fr["sstype"] == "ah" and len(fr["reslist"]) >= min_ah) or (
                              fr["sstype"] == "bs" and len(fr["reslist"]) >= min_bs)]


            add_coil = []
            listaFrags = sorted(listaFrags,key=lambda x: x["reslist"][-1])
            for o, frag in enumerate(listaFrags):
                if o > 0 and len(listaFrags[o - 1]["reslist"]) > 0 and len(listaFrags[o]["reslist"]) > 0 and \
                                listaFrags[o - 1]["sstype"] in ["ah", "bs"] and listaFrags[o]["sstype"] in ["ah", "bs"]:
                    resaN = Bioinformatics.get_residue(strucc, frag["reslist"][0][1], frag["reslist"][0][2], frag["reslist"][0][3])
                    prev = listaFrags[o - 1]["reslist"][-1]
                    prevResC = Bioinformatics.get_residue(strucc, prev[1], prev[2], prev[3])
                    #print(o, prev, frag["reslist"][0], Bioinformatics.check_continuity(resaN, prevResC))
                    #geome = compute_instruction(get_unique_cv_among_residues(strucc, listaFrags[o - 1]["reslist"]), get_unique_cv_among_residues(strucc, listaFrags[o]["reslist"]), unique_fragment_cv=True)

                    if Bioinformatics.check_continuity(resaN, prevResC):
                        #print([resi[4] for resi in listaFrags[o - 1]["reslist"][-3:]])
                        #print([resi[4] for resi in listaFrags[o]["reslist"][:3]])
                        try:
                            id1_gly = [resi[4] for resi in listaFrags[o - 1]["reslist"][-3:]].index("GLY")
                        except:
                            id1_gly = -1
                        try:
                            id1_pro = [resi[4] for resi in listaFrags[o - 1]["reslist"][-3:]].index("PRO")
                        except:
                            id1_pro = -1
                        try:
                            id2_gly = [resi[4] for resi in listaFrags[o]["reslist"][:3]].index("GLY")
                        except:
                            id2_gly = -1
                        try:
                            id2_pro = [resi[4] for resi in listaFrags[o]["reslist"][:3]].index("PRO")
                        except:
                            id2_pro = -1
                        #print(id1_gly,id1_pro,id2_gly,id2_pro)

                        A = False
                        B = False
                        if id1_gly >= 0:
                            dirco = {"sstype": "coil", "reslist": [listaFrags[o - 1]["reslist"][-3+id1_gly]], "cvids": [listaFrags[o - 1]["cvids"][-1]] ,
                                     "cvls": [listaFrags[o - 1]["cvls"][-1]], "sequence": "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], [listaFrags[o - 1]["reslist"][-3 + id1_gly]])))}

                            add_coil.append(dirco)

                            if id1_gly < 2:
                                listaFrags[o]["reslist"] = listaFrags[o - 1]["reslist"][-3+id1_gly+1:] + listaFrags[o]["reslist"]
                                listaFrags[o]["sequence"] =  "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))
                                #listaFrags[o]["cvids"] = [listaFrags[o - 1]["cvids"][-1]] + listaFrags[o]["cvids"]
                                #listaFrags[o]["cvls"] = [listaFrags[o - 1]["cvls"][-1]] + listaFrags[o]["cvls"]
                            listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"][:-3+id1_gly]
                            listaFrags[o - 1]["sequence"] = "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))
                            #listaFrags[o - 1]["cvids"] =  listaFrags[o - 1]["cvids"][:-1]
                            #listaFrags[o - 1]["cvls"] = listaFrags[o - 1]["cvls"][:-1]

                            A = True
                            #print("A",listaFrags[o - 1]["reslist"])
                            #print("A",listaFrags[o]["reslist"])
                        elif id1_pro >= 0:
                            dirco = {"sstype": "coil", "reslist": [listaFrags[o - 1]["reslist"][-3 + id1_pro]],
                                     "cvids": [listaFrags[o - 1]["cvids"][-1]],
                                     "cvls": [listaFrags[o - 1]["cvls"][-1]], "sequence": "".join(list(
                                    map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                        [listaFrags[o - 1]["reslist"][-3 + id1_pro]])))}

                            add_coil.append(dirco)

                            if id1_pro < 2:
                                listaFrags[o]["reslist"] = listaFrags[o - 1]["reslist"][-3+id1_pro+1:] + listaFrags[o]["reslist"]
                                listaFrags[o]["sequence"] =  "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))
                            listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"][:-3+id1_pro]
                            listaFrags[o - 1]["sequence"] = "".join(
                                list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                            A = True
                            #print("B", listaFrags[o - 1]["reslist"])
                            #print("B", listaFrags[o]["reslist"])
                        elif id2_gly >= 0:
                            dirco = {"sstype": "coil", "reslist": [listaFrags[o]["reslist"][id2_gly]],
                                     "cvids": [listaFrags[o]["cvids"][-1]],
                                     "cvls": [listaFrags[o]["cvls"][-1]], "sequence": "".join(list(
                                    map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                        [listaFrags[o]["reslist"][id2_gly]])))}

                            add_coil.append(dirco)
                            listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"] + listaFrags[o]["reslist"][:id2_gly]
                            listaFrags[o - 1]["sequence"] = "".join(
                                list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                            if len(listaFrags[o]["reslist"]) > 3:
                                listaFrags[o]["reslist"] = listaFrags[o]["reslist"][id2_gly+1:]
                                listaFrags[o]["sequence"] =  "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))

                            B = True
                            #print("C", listaFrags[o - 1]["reslist"])
                            #print("C", listaFrags[o]["reslist"])
                        elif id2_pro >= 0:
                            dirco = {"sstype": "coil", "reslist": [listaFrags[o]["reslist"][id2_pro]],
                                     "cvids": [listaFrags[o]["cvids"][-1]],
                                     "cvls": [listaFrags[o]["cvls"][-1]], "sequence": "".join(list(
                                    map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                        [listaFrags[o]["reslist"][id2_pro]])))}

                            add_coil.append(dirco)
                            listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"] + listaFrags[o]["reslist"][:id2_pro]
                            listaFrags[o - 1]["sequence"] = "".join(
                                list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                            if len(listaFrags[o]["reslist"]) > 3:
                                listaFrags[o]["reslist"] = listaFrags[o]["reslist"][id2_pro + 1:]
                                listaFrags[o]["sequence"] =  "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))

                            B = True
                            #print("D", listaFrags[o - 1]["reslist"])
                            #print("D", listaFrags[o]["reslist"])
                        else:
                            dirco = {"sstype": "coil", "reslist": [listaFrags[o - 1]["reslist"][-1]],
                                     "cvids": [listaFrags[o - 1]["cvids"][-1]],
                                     "cvls": [listaFrags[o - 1]["cvls"][-1]], "sequence": "".join(list(
                                    map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                        [listaFrags[o - 1]["reslist"][-1]])))}

                            add_coil.append(dirco)
                            listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"][:-1]
                            listaFrags[o - 1]["sequence"] = "".join(
                                list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                            A = True
                            #print("E", listaFrags[o - 1]["reslist"])
                            #print("E", listaFrags[o]["reslist"])

                        #print(add_coil)
                        #quit()

                        if A and len(listaFrags[o-1]["reslist"]) < 3:
                            listaFrags[o - 1]["sstype"] = "coil"
                        elif B and len(listaFrags[o]["reslist"]) < 3:
                            listaFrags[o]["sstype"] = "coil"

            listaFrags = [fr for fr in listaFrags if len(fr["reslist"])>0]
            listaFrags = sorted(listaFrags,key=lambda x: x["reslist"][-1])

        g = igraph.Graph.Full(len(listaFrags))

        for i, fr in enumerate(listaFrags):
            g.vs[i]["sstype"] = fr["sstype"]
            g.vs[i]["reslist"] = fr["reslist"]
            g.vs[i]["unique_cv"] = get_unique_cv_among_residues(strucc, fr["reslist"])
            try:
                g.vs[i]["vecLength"] = get_atoms_distance(g.vs[i]["unique_cv"][2],g.vs[i]["unique_cv"][3])
            except:
                g.vs[i]["vecLength"] = -1
            g.vs[i]["cvids"] = fr["cvids"]
            g.vs[i]["resforgraph"] = fr["resforgraph"]
            g.vs[i]["seqforgraph"] = fr["seqforgraph"]
            # indfri = int(len(fr["cvids"]) / 2)
            # indfr = fr["cvids"][indfri]
            # cvst = list(filter(lambda x: x[0] == indfr, cvs_list))[0]
            # g.vs[i]["com"] = [(cvst[2][0] + cvst[3][0]) / 2.0, (cvst[2][1] + cvst[3][1]) / 2.0,
            #                   (cvst[2][2] + cvst[3][2]) / 2.0]
            g.vs[i]["cvls"] = fr["cvls"]
            g.vs[i]["sequence"] = fr["sequence"]
        return g, a

    def __generate_graph(a, min_ah, min_bs, associations, interface = False):
        if min_ah is not None and min_ah < 3:
            min_ah = 3
        if min_bs is not None and min_bs < 3:
            min_bs = 3

        listaFrags = []

        if interface:
            associations_sorted = sorted(associations.keys(), key=lambda x: x[:3]+(x[3][2],)+(x[3][1],)) #Sorted by chain, then insertion_code and then residue
        else:
            associations_sorted = sorted(associations.keys(), key=lambda x: x[:3]+(x[3][1:],))

        # for asso in sorted(associations.keys(), key=lambda x: x[:3]+(x[3][1:],)): #The way it was before interfaces
        for asso in associations_sorted:
            if len(listaFrags) == 0:
                listaFrags.append([])
            if len(listaFrags[-1]) == 0:
                listaFrags[-1].append(asso)
            elif len(listaFrags[-1]) > 0 and associations[listaFrags[-1][-1]]["result"].lower() == associations[asso]["result"].lower():
                resi = asso
                prev = listaFrags[-1][-1]
                resaN = Bioinformatics.get_residue(strucc, resi[1], resi[2], resi[3])
                prevResC = Bioinformatics.get_residue(strucc, prev[1], prev[2], prev[3])
                if Bioinformatics.check_continuity(resaN, prevResC):
                    listaFrags[-1].append(asso)
                else:
                    listaFrags.append([asso])
            else:
                listaFrags.append([asso])

        listaFrags = [
            {"sstype": associations[fragl[0]]["result"].lower(),
             "reslist": [list(fr) for fr in fragl],
             "cvids": sorted(set([e[3] for fr in fragl for e in a if list(fr) in e[2]])),
             "sequence": "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], fragl))),
             "strictnesses": [t[1] for t in sorted(set([(e[3],e[4]) for fr in fragl for e in a if list(fr) in e[2] if len(e)==5]), key=lambda hl: hl[0])]}
            for fragl in listaFrags if len(fragl) > 0 and associations[fragl[0]]["result"].lower() in ["ah", "bs", "coil", "COIL"]]

        for o,frag in enumerate(listaFrags):
            if o > 0 and len(listaFrags[o - 1]["reslist"]) > 0 and len(set(listaFrags[o - 1]["cvids"])&set(frag["cvids"])) > 0:
                listaFrags[o - 1]["cvids"] = sorted(list(set(listaFrags[o - 1]["cvids"])-set(frag["cvids"])))
                listaFrags[o - 1]["cvls"] = [e[0] for e in a if e[3] in listaFrags[o - 1]["cvids"]]
                # print(set([(ss,listaFrags[o - 1]["strictnesses"][pp]) for pp,ss in enumerate(listaFrags[o - 1]["cvids"])]))
                # print(set([(ss,frag["strictnesses"][pp]) for pp,ss in enumerate(frag["cvids"])]))
                # print(list(set([(ss,listaFrags[o - 1]["strictnesses"][pp]) for pp,ss in enumerate(listaFrags[o - 1]["cvids"])])-set([(ss,frag["strictnesses"][pp]) for pp,ss in frag["cvids"]])))
                # print(sorted(list(set([(ss,listaFrags[o - 1]["strictnesses"][pp]) for pp,ss in enumerate(listaFrags[o - 1]["cvids"])])-set([(ss,frag["strictnesses"][pp]) for pp,ss in frag["cvids"]])), key=lambda hl: hl[0]))
                # print([t[1] for t in sorted(list(set([(ss,listaFrags[o - 1]["strictnesses"][pp]) for pp,ss in enumerate(listaFrags[o - 1]["cvids"])])-set([(ss,frag["strictnesses"][pp]) for pp,ss in frag["cvids"]])), key=lambda hl: hl[0])])
                listaFrags[o - 1]["strictnesses"] = [t[1] for t in sorted(list(set([(ss,listaFrags[o - 1]["strictnesses"][pp]) for pp,ss in enumerate(listaFrags[o - 1]["cvids"])])-set([(ss,frag["strictnesses"][pp]) for pp,ss in enumerate(frag["cvids"])])), key=lambda hl: hl[0])]
            listaFrags[o]["cvls"] = [e[0] for e in a if e[3] in listaFrags[o]["cvids"]]

        if min_ah is not None and min_bs is not None:
            # print('min_ah,min_bs',min_ah,min_bs)

            listaFrags = [fr for i, fr in enumerate(listaFrags) if
                          (fr["sstype"] == "coil" and len(fr["reslist"]) > 0) or
                          (fr["sstype"] == "ah" and len(fr["reslist"]) >= min_ah) or (
                                  fr["sstype"] == "bs" and len(fr["reslist"]) >= min_bs)]
            # for frag in listaFrags:
            #     if frag['sstype']!='coil':
            #         print("frag['sstype']",frag['sstype'],"len(frag['reslist']",len(frag['reslist']))

        add_coil = []
        listaFrags = sorted(listaFrags, key=lambda x: x["reslist"][-1])
        for o, frag in enumerate(listaFrags):
            if o > 0 and len(listaFrags[o - 1]["reslist"]) > 2 and len(listaFrags[o]["reslist"]) > 2 and \
                    listaFrags[o - 1]["sstype"] in ["ah", "bs"] and listaFrags[o]["sstype"] in ["ah", "bs"]:
                resaN = Bioinformatics.get_residue(strucc, frag["reslist"][0][1], frag["reslist"][0][2], frag["reslist"][0][3])
                prev = listaFrags[o - 1]["reslist"][-1]
                prevResC = Bioinformatics.get_residue(strucc, prev[1], prev[2], prev[3])

                if Bioinformatics.check_continuity(resaN, prevResC):
                    # print([resi[4] for resi in listaFrags[o - 1]["reslist"][-3:]])
                    # print([resi[4] for resi in listaFrags[o]["reslist"][:3]])
                    try:
                        id1_gly = [resi[4] for resi in listaFrags[o - 1]["reslist"][-3:]].index("GLY")
                    except:
                        id1_gly = -1
                    try:
                        id1_pro = [resi[4] for resi in listaFrags[o - 1]["reslist"][-3:]].index("PRO")
                    except:
                        id1_pro = -1
                    try:
                        id2_gly = [resi[4] for resi in listaFrags[o]["reslist"][:3]].index("GLY")
                    except:
                        id2_gly = -1
                    try:
                        id2_pro = [resi[4] for resi in listaFrags[o]["reslist"][:3]].index("PRO")
                    except:
                        id2_pro = -1
                    # print(id1_gly,id1_pro,id2_gly,id2_pro)

                    A = False
                    B = False
                    if id1_gly >= 0:
                        dirco = {"sstype": "coil", "reslist": [listaFrags[o - 1]["reslist"][-3 + id1_gly]],
                                 "cvids": [listaFrags[o - 1]["cvids"][-1]],
                                 "cvls": [listaFrags[o - 1]["cvls"][-1]], "sequence": "".join(list(
                                map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                    [listaFrags[o - 1]["reslist"][-3 + id1_gly]]))),
                                 "strictnesses": [listaFrags[o - 1]["strictnesses"][-1]]}

                        add_coil.append(dirco)

                        if id1_gly < 2:
                            listaFrags[o]["reslist"] = listaFrags[o - 1]["reslist"][-3 + id1_gly + 1:] + listaFrags[o]["reslist"]
                            listaFrags[o]["sequence"] = "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))
                            # listaFrags[o]["cvids"] = [listaFrags[o - 1]["cvids"][-1]] + listaFrags[o]["cvids"]
                            # listaFrags[o]["cvls"] = [listaFrags[o - 1]["cvls"][-1]] + listaFrags[o]["cvls"]
                        listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"][:-3 + id1_gly]
                        listaFrags[o - 1]["sequence"] = "".join(
                            list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))
                        # listaFrags[o - 1]["cvids"] =  listaFrags[o - 1]["cvids"][:-1]
                        # listaFrags[o - 1]["cvls"] = listaFrags[o - 1]["cvls"][:-1]

                        A = True
                        # print("A",listaFrags[o - 1]["reslist"])
                        # print("A",listaFrags[o]["reslist"])
                    elif id1_pro >= 0:
                        dirco = {"sstype": "coil", "reslist": [listaFrags[o - 1]["reslist"][-3 + id1_pro]],
                                 "cvids": [listaFrags[o - 1]["cvids"][-1]],
                                 "cvls": [listaFrags[o - 1]["cvls"][-1]], "sequence": "".join(list(
                                map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                    [listaFrags[o - 1]["reslist"][-3 + id1_pro]]))),
                                 "strictnesses": [listaFrags[o - 1]["strictnesses"][-1]]}

                        add_coil.append(dirco)

                        if id1_pro < 2:
                            listaFrags[o]["reslist"] = listaFrags[o - 1]["reslist"][-3 + id1_pro + 1:] + listaFrags[o]["reslist"]
                            listaFrags[o]["sequence"] = "".join(
                                list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))
                        listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"][:-3 + id1_pro]
                        listaFrags[o - 1]["sequence"] = "".join(
                            list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                        A = True
                        # print("B", listaFrags[o - 1]["reslist"])
                        # print("B", listaFrags[o]["reslist"])
                    elif id2_gly >= 0:
                        dirco = {"sstype": "coil", "reslist": [listaFrags[o]["reslist"][id2_gly]],
                                 "cvids": [listaFrags[o]["cvids"][-1]],
                                 "cvls": [listaFrags[o]["cvls"][-1]], "sequence": "".join(list(
                                map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                    [listaFrags[o]["reslist"][id2_gly]]))),
                                 "strictnesses": [listaFrags[o]["strictnesses"][-1]]}

                        add_coil.append(dirco)
                        listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"] + listaFrags[o]["reslist"][
                                                                                      :id2_gly]
                        listaFrags[o - 1]["sequence"] = "".join(
                            list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                        if len(listaFrags[o]["reslist"]) > 3:
                            listaFrags[o]["reslist"] = listaFrags[o]["reslist"][id2_gly + 1:]
                            listaFrags[o]["sequence"] = "".join(
                                list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))

                        B = True
                        # print("C", listaFrags[o - 1]["reslist"])
                        # print("C", listaFrags[o]["reslist"])
                    elif id2_pro >= 0:
                        dirco = {"sstype": "coil", "reslist": [listaFrags[o]["reslist"][id2_pro]],
                                 "cvids": [listaFrags[o]["cvids"][-1]],
                                 "cvls": [listaFrags[o]["cvls"][-1]], "sequence": "".join(list(
                                map(lambda x: Bioinformatics.AADICMAP[x[4]],
                                    [listaFrags[o]["reslist"][id2_pro]]))),
                                 "strictnesses": [listaFrags[o]["strictnesses"][-1]]}

                        add_coil.append(dirco)
                        listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"] + listaFrags[o]["reslist"][:id2_pro]
                        listaFrags[o - 1]["sequence"] = "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                        if len(listaFrags[o]["reslist"]) > 3:
                            listaFrags[o]["reslist"] = listaFrags[o]["reslist"][id2_pro + 1:]
                            listaFrags[o]["sequence"] = "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o]["reslist"])))

                        B = True
                        # print("D", listaFrags[o - 1]["reslist"])
                        # print("D", listaFrags[o]["reslist"])
                    else:
                        dirco = {"sstype": "coil", "reslist": [listaFrags[o - 1]["reslist"][-1]],
                                 "cvids": [listaFrags[o - 1]["cvids"][-1]],
                                 "cvls": [listaFrags[o - 1]["cvls"][-1]], "sequence": "".join(list(map(lambda x: Bioinformatics.AADICMAP[x[4]], [listaFrags[o - 1]["reslist"][-1]]))),
                                 "strictnesses": [listaFrags[o - 1]["strictnesses"][-1]]}

                        add_coil.append(dirco)
                        listaFrags[o - 1]["reslist"] = listaFrags[o - 1]["reslist"][:-1]
                        listaFrags[o - 1]["sequence"] = "".join(
                            list(map(lambda x: Bioinformatics.AADICMAP[x[4]], listaFrags[o - 1]["reslist"])))

                        A = True
                        # print("E", listaFrags[o - 1]["reslist"])
                        # print("E", listaFrags[o]["reslist"])

                    # print(add_coil)
                    # quit()

                    if A and len(listaFrags[o - 1]["reslist"]) < 3:
                        listaFrags[o - 1]["sstype"] = "coil"
                    elif B and len(listaFrags[o]["reslist"]) < 3:
                        listaFrags[o]["sstype"] = "coil"

        listaFrags = [fr for fr in listaFrags if len(fr["reslist"]) > 0]
        listaFrags = sorted(listaFrags, key=lambda x: x["reslist"][-1])

        g = igraph.Graph.Full(len(listaFrags))

        for i, fr in enumerate(listaFrags):
            g.vs[i]["sstype"] = fr["sstype"]
            g.vs[i]["reslist"] = fr["reslist"]
            g.vs[i]["unique_cv"] = get_unique_cv_among_residues(strucc, fr["reslist"])
            try:
                g.vs[i]["vecLength"] = get_atoms_distance(g.vs[i]["unique_cv"][2], g.vs[i]["unique_cv"][3])
            except:
                g.vs[i]["vecLength"] = -1
            g.vs[i]["cvids"] = fr["cvids"]
            #g.vs[i]["resforgraph"] = fr["reslist"]
            #g.vs[i]["seqforgraph"] = fr["sequence"]
            # indfri = int(len(fr["cvids"]) / 2)
            # indfr = fr["cvids"][indfri]
            # cvst = list(filter(lambda x: x[0] == indfr, cvs_list))[0]
            # g.vs[i]["com"] = [(cvst[2][0] + cvst[3][0]) / 2.0, (cvst[2][1] + cvst[3][1]) / 2.0,
            #                   (cvst[2][2] + cvst[3][2]) / 2.0]
            g.vs[i]["cvls"] = fr["cvls"]
            g.vs[i]["sequence"] = fr["sequence"]
            g.vs[i]["strictnesses"] = fr["strictnesses"]

        # print("===============================LISTAFRAGS===============================")
        # for fra in listaFrags:
        #      print(fra["sstype"],fra['sequence'])
        # print("=========================================================================")
        return g, a

    def __generate_3d_relations(g, max_distance=10.0, validate=["bs", "coil"]):
        dictio3d = {}
        for i, frag1 in enumerate(g.vs):
            if frag1["sstype"] not in validate:
                continue
            for frag2 in g.vs.select(lambda vertex: vertex.index != frag1.index):
                if frag2["sstype"] not in validate:
                    continue
                for uno in frag1["cvids"]:
                    if uno not in dictio3d:
                        dictio3d[uno] = []
                    for due in frag2["cvids"]:
                        if due == uno:
                            continue
                        # print "||", uno, "--", due,"||",
                        ###value1 = compute_instruction(cvs_list[dimap[uno]], cvs_list[dimap[due]])
                        sup = tuple(sorted([dimap[uno], dimap[due]]))
                        value1 = matrix[sup]
                        if value1[3] <= max_distance:
                            dictio3d[uno].append([due, value1[2], value1[3], value1[4]])
        return dictio3d

    def __check_impossible_angle(z,p,a,enne,zenne,cvs_list,dimap,matrix,angle_mean_bs,angle_mean_ah,strucc):
        # print("========CHECKING=========")
        # if p - 3 >= 0:
        #     print(enne)
        #     print(cvs_list[dimap[z[p - 3][3]]])
        #     print(abs(matrix[tuple(sorted([dimap[enne[3]], dimap[z[p - 3][3]]]))][2] - angle_mean_bs))
        # print()
        # if p - 2 >= 0:
        #     print(zenne)
        #     print(cvs_list[dimap[z[p - 2][3]]])
        #     print(abs(matrix[tuple(sorted([dimap[zenne[3]], dimap[z[p - 2][3]]]))][2] - angle_mean_bs))
        # print("==========================")
        # Checking if the current enne is annotated as ah or bs and if it exists a previous cv at p-3 with the same annotation
        if enne[1] in ["ah", "bs"] and p - 3 >= 0 and z[p - 3][1] == enne[1]:
            # checking if the two CV are representing two tripetides continous but not overlapping
            resaN = Bioinformatics.get_residue(strucc, enne[2][0][1], enne[2][0][2], enne[2][0][3])
            prevResC = Bioinformatics.get_residue(strucc, z[p - 3][2][-1][1], z[p - 3][2][-1][2], z[p - 3][2][-1][3])
            if Bioinformatics.check_continuity(resaN, prevResC):
                # checking if their angle is not too far away from the standard mean of the annotation
                if (enne[1] == "bs" and abs(
                        matrix[tuple(sorted([dimap[enne[3]], dimap[z[p - 3][3]]]))][2] - angle_mean_bs) > 50) \
                        or (enne[1] == "ah" and abs(
                    matrix[tuple(sorted([dimap[enne[3]], dimap[z[p - 3][3]]]))][2] - angle_mean_ah) > 50):
                    # print("We found a case where enne is")
                    # print(enne)
                    # print("And z[p-3] is")
                    # print(z[p-3])
                    # print("the angle between them is",matrix[tuple(sorted([dimap[enne[3]], dimap[z[p-3][3]]]))][2])
                    # enne[1] = "COIL"
                    if z[p - 3][1] not in ["coil", "COIL"]:
                        z[p - 2][1] = "COIL"
                    else:
                        z[p - 3][1] = "COIL"
        # Checking if the current zenne is annotated as ah or bs and if it exists a previous cv at p-2 with the same annotation
        if zenne[1] in ["ah", "bs"] and p - 2 >= 0 and z[p - 2][1] == zenne[1]:
            # checking if the two CV are representing two tripetides continous but not overlapping
            resaN = Bioinformatics.get_residue(strucc, zenne[2][0][1], zenne[2][0][2], zenne[2][0][3])
            prevResC = Bioinformatics.get_residue(strucc, z[p - 2][2][-1][1], z[p - 2][2][-1][2], z[p - 2][2][-1][3])
            if Bioinformatics.check_continuity(resaN, prevResC):
                # checking if their angle is not too far away from the standard mean of the annotation
                if (zenne[1] == "bs" and abs(
                        matrix[tuple(sorted([dimap[zenne[3]], dimap[z[p - 2][3]]]))][2] - angle_mean_bs) > 50) \
                        or (zenne[1] == "ah" and abs(
                    matrix[tuple(sorted([dimap[zenne[3]], dimap[z[p - 2][3]]]))][2] - angle_mean_ah) > 50):
                    # print("We found a case where zenne is")
                    # print(zenne)
                    # print("And z[p-2] is")
                    # print(z[p - 2])
                    # print("the angle between them is", matrix[tuple(sorted([dimap[zenne[3]], dimap[z[p - 2][3]]]))][2])
                    # zenne[1] = "COIL"
                    if z[p - 2][1] not in ["coil", "COIL"]:
                        z[p - 1][1] = "COIL"
                    else:
                        z[p - 2][1] = "COIL"
        z[p] = enne
        z[p + 1] = zenne
        return z

    z = [0 for y in range(len(a))]
    tf = True
    for p in range(len(a) - 1):
        enne, zenne, tf = __check_by_unified_score_step2(a[p], a[p + 1], None, take_first=tf)
        z = __check_impossible_angle(z, p, a, enne, zenne, cvs_list, dimap, matrix, angle_mean_bs, angle_mean_ah, strucc)
    a = z

    text = ""
    for w in range(2):
        associations = __get_associations(a, stringent=True if w == 0 else False)
        # associations = __getAssociations(a)
        g, a = __generate_graph(a, None if w == 0 else min_ah, None if w == 0 else min_bs, associations, interface)
        # for frag in g.vs:
        #      if len(frag["reslist"]) > 0:
        #          print(frag["reslist"][0][3][1], "--", frag["reslist"][-1][3][1], frag["sstype"], "---", frag["unique_cv"])

        dictio_3D = __generate_3d_relations(g, validate=["bs", "coil"] if w == 0 else ["bs"])

        #NOTE: Why I was filtering here by the length of e the problem.
        #a = [e[:-1] if len(e) == 5 else e for e in a]

        z = [0 for y in range(len(a))]
        tf = True
        for p in range(len(a) - 1):
            enne, zenne, tf = __check_by_unified_score_step2(a[p], a[p + 1], dictio_3D, take_first=tf,
                                                             validate=["bs", "coil"] if w == 0 else ["bs"],
                                                             min_num_bs=1.0 if w == 0 else 1.0)

            z = __check_impossible_angle(z, p, a, enne, zenne, cvs_list, dimap, matrix, angle_mean_bs, angle_mean_ah, strucc)

        a = z

        # for e in a:
        #     print(e)
        # print()
        # print()

        if w == 1:
            associations = __get_associations(a, stringent=True if w == 0 else False)
            # associations = __getAssociations(a)
            g, a = __generate_graph(a, min_ah, min_bs, associations, interface = interface)

        for e in a:
            text += str(e)+"\n"
        text += "\n"
    # for frag in listaFrags:
    #     print([[(w[2],w[3][1]) for w in fr[2]] for fr in frag],len([fr[2] for fr in frag]),[fr[3] for fr in frag],len([fr[3] for fr in frag]))
    #     print([[(w[2],w[3][1]) for w in fr[2]] for fr in frag[:1]]+[[(w[2],w[3][1]) for w in fr[2][-1:]] for fr in frag[1:]],len([[(w[2],w[3][1]) for w in fr[2]] for fr in frag[:1]]+[[(w[2],w[3][1]) for w in fr[2][-2:-1]] for fr in frag[1:]]),[fr[3] for fr in frag],len([fr[3] for fr in frag]))
    #     print(frag)
    #     print()

    with open("fromCVtoAA_12.txt","w") as f:
        f.write(text)

    return g

@SystemUtility.timing
# @SystemUtility.profileit
def aleph_internal_terstr(strucc, g, matrix, cvs):
    """

    :param strucc:
    :type strucc:
    :param g:
    :type g:
    :param cvs:
    :type cvs:
    :return:
    :rtype:
    """
    dimap = {value[0]: i for (i, value) in enumerate(cvs)}

    for i, frag1 in enumerate(g.vs):
        if len(frag1["cvids"]) > 0:
            alls = numpy.empty((len(frag1["cvids"]) - 1, 3))
            for u in range(len(frag1["cvids"]) - 1):
                a = frag1["cvids"][u]
                b = frag1["cvids"][u + 1]
                ###value1 = compute_instruction(cvs[dimap[a]], cvs[dimap[b]])
                sup = tuple(sorted([dimap[a], dimap[b]]))
                value1 = matrix[sup]
                alls[u] = numpy.array([value1[2], value1[3], value1[4]])
            frag1["alls"] = alls
    return g

@SystemUtility.timing
# @SystemUtility.profileit
def aleph_terstr(strucc, g, matrix, cvs, weight="distance_avg", verbose=True):
    """

    :param strucc:
    :type strucc:
    :param g:
    :type g:
    :param cvs:
    :type cvs:
    :param weight:
    :type weight:
    :param full:
    :type full:
    :param sheet_factor:
    :type sheet_factor:
    :return:
    :rtype:
    """

    # print "LEN_CVS:",len(cvs)
    dimap = {value[0]: i for (i, value) in enumerate(cvs)}
    sstype_map = {"ah": 1, "bs": 2, "coil": 3}

    convert = {}
    for i, frag1 in enumerate(g.vs):
        framme1 = []
        for frag2 in g.vs.select(lambda vertex: vertex.index > frag1.index):
            # g.add_edge(frag1.index, frag2.index)
            t = 0.0
            z = [1.0, 1.0, 1.0]
            minimum = 10000000000000.0
            alls = []
            minforcvids = []
            #maximum_number_of_mins = len(frag1["cvids"]) #min(len(frag1["cvids"]),len(frag2["cvids"]))
            fgs = sorted([frag1,frag2],key=lambda x: len(x["cvids"]))
            for a in fgs[0]["cvids"]:
                minoide = [[0.0,100.0,0]]
                for b in fgs[1]["cvids"]:
                    # print "test tertiary frags:",i,"====","a,b,c",a,b,c,"d,e,f",d,e,f
                    ###value1 = compute_instruction(cvs[dimap[a]], cvs[dimap[b]])
                    sup = tuple(sorted([dimap[a], dimap[b]]))
                    value1 = matrix[sup]
                    z = [z[0] + value1[2], z[1] + value1[3], z[2] + value1[4]]
                    minimum = min(minimum, value1[3])
                    alls.append([a, b, value1[2], value1[3], value1[4]])
                    minoide.append([value1[2], value1[3], value1[4]])
                    t += 1
                minforcvids.append(sorted(minoide, key=lambda x: x[1])[0])

            if t==0: t=1
            z = [z[0] / t, z[1] / t, z[2] / t, sstype_map[frag1["sstype"]], sstype_map[frag2["sstype"]]]
            #print(frag1)
            #print(frag2)
            g.es[g.get_eid(frag1.index, frag2.index)]["mean"] = z
            if frag1["sstype"] in ["ah","bs"] and frag2["sstype"] in ["ah","bs"]:
                g.es[g.get_eid(frag1.index, frag2.index)]["mean_cv_uniques"] = compute_instruction(g.vs[frag1.index]["unique_cv"],g.vs[frag2.index]["unique_cv"],unique_fragment_cv=True)[2:6]+z[-2:]
            g.es[g.get_eid(frag1.index, frag2.index)]["avg"] = z[1]
            g.es[g.get_eid(frag1.index, frag2.index)]["min"] = minimum
            g.es[g.get_eid(frag1.index, frag2.index)]["alls"] = alls

            if weight == "distance_avg":
                if frag1["sstype"] == "bs" and frag2["sstype"] == "bs":
                    g.es[g.get_eid(frag1.index, frag2.index)]["weight"] = 100*(100.0/z[1])
                    g.es[g.get_eid(frag1.index, frag2.index)]["spantree_weight"] = z[1]/100.0
                elif frag1["sstype"] == "ah" and frag2["sstype"] == "ah":
                    g.es[g.get_eid(frag1.index, frag2.index)]["weight"] = 80 * (100.0 / z[1])
                    g.es[g.get_eid(frag1.index, frag2.index)]["spantree_weight"] = z[1] / 80.0
                else:
                    g.es[g.get_eid(frag1.index, frag2.index)]["weight"] =  (100.0 / z[1])
                    g.es[g.get_eid(frag1.index, frag2.index)]["spantree_weight"] = z[1]
            elif weight == "distance_min":
                if frag1["sstype"] == "bs" and frag2["sstype"] == "bs":
                    g.es[g.get_eid(frag1.index, frag2.index)]["weight"] = 100 * (100.0 / minimum)
                    g.es[g.get_eid(frag1.index, frag2.index)]["spantree_weight"] = minimum / 100.0
                elif frag1["sstype"] == "ah" and frag2["sstype"] == "ah":
                    g.es[g.get_eid(frag1.index, frag2.index)]["weight"] = 80 * (100.0 / minimum)
                    g.es[g.get_eid(frag1.index, frag2.index)]["spantree_weight"] = minimum / 80.0
                else:
                    g.es[g.get_eid(frag1.index, frag2.index)]["weight"] = (100.0 / minimum)
                    g.es[g.get_eid(frag1.index, frag2.index)]["spantree_weight"] = minimum

            if frag1["sstype"] == "bs" and frag2["sstype"] == "bs" and len(minforcvids) > 0:
                teren = len(minforcvids) if len(minforcvids)<3 else len(minforcvids)-2
                q = (sum([max(BS_UD_EA[int(round(ed[0]))],BS_UU_EA[int(round(ed[0]))])/BS_MAX[numpy.argmax([BS_UD_EA[int(round(ed[0]))],BS_UU_EA[int(round(ed[0]))]])] for ed in minforcvids if ed[1] <= 6])/(teren))*100.0
                framme1.append((q,frag2))

        if frag1["sstype"] == "bs":
            for q,frag2 in sorted(framme1,key=lambda x: x[0], reverse=True)[:2]:
                if verbose: print("PERCENTAGE OF SHEET", q, frag1["sequence"], frag1["reslist"][0][2], frag2["sequence"], frag2["reslist"][0][2], len(frag1["cvids"]), len(frag2["cvids"]))
                if q >= 35:
                    if frag1["sheet"] is not None:
                        if frag2["sheet"] is not None:
                            convert[frag2["sheet"]] = frag1["sheet"]
                        frag2["sheet"] = frag1["sheet"]
                    elif frag2["sheet"] is not None:
                        if frag1["sheet"] is not None:
                             convert[frag1["sheet"]] = frag2["sheet"]
                        frag1["sheet"] = frag2["sheet"]
                    else:
                        frag1["sheet"] = max([u for u in g.vs["sheet"] if u is not None])+1 if len([u for u in g.vs["sheet"] if u is not None]) > 0 else 0
                        frag2["sheet"] = frag1["sheet"]
                        #print("---",frag1["sequence"],frag2["sequence"],frag1["sheet"],frag2["sheet"])
                # else:
                #     if frag1["sheet"] is not None and frag2["sheet"] is not None and frag1["sheet"] == frag2["sheet"]:
                #         frag2["sheet"] = None

    for key in sorted(convert.keys()):
        g.vs["sheet"] = [n if n is None else n if n != key else convert[key]  for n in  g.vs["sheet"]]

    if verbose:
        for frag in g.vs:
           if frag["sstype"] == "bs" and frag["sheet"] is None:
               #print("The fragment",frag["sstype"],frag["sequence"],"is converted to coil because is not packed in a sheet")
               #frag["sstype"] = "coil"
               print("The fragment", frag["sstype"], frag["sequence"],"has None sheet id","chain",frag["reslist"][0][2])
           elif frag["sstype"] == "bs":
               print("The fragment",frag["sstype"],frag["sequence"],"has sheet id:",frag["sheet"],"chain",frag["reslist"][0][2])

    return g

def annotate_pdb_model_with_aleph(pdb_model, structure=None, weight="distance_avg", min_ah=4, min_bs=3, write_pdb=True,
                                  strictness_ah=0.45, strictness_bs=0.20, peptide_length=3, is_model=False, only_reformat=True,
                                  path_base="./", verbose=True, process_only_chains=None, interface=False, domains=False,
                                  compactness=False):
    """
     Annotates a protein pdb file with CVs and builds a graph in which each node is a secondary structure element
    or a coil fragment and the edge connecting two nodes is the metric distance between these fragments.

    :param pdb_model: The pdb file to annotate
    :type pdb_model: io.TextIOWrapper
    :param weight: The weight scheme for computing edge parameters. ["distance_avg","distance_min","distance_max"]
    :type weight: str
    :param min_ah: Minimum number of residues for an alpha helix to be accepted
    :type min_ah: int
    :param min_bs: Minimum number of residues for a beta strand to be accepted
    :type min_bs: int
    :param write_pdb: Write an annotated pdb file
    :type write_pdb: bool
    :param strictness_ah: strictness parameter threshold for accepting ah CVs
    :type strictness_ah: float
    :param strictness_bs: strictness parameter threshold for accepting bs CVs
    :type strictness_bs: float
    :param peptide_length: Define the peptide length for computing a CV
    :type peptide_length: int
    :param is_model:
    :type is_model:
    :param only_reformat:
    :type only_reformat:
    :param path_base:
    :type path_base:
    :return
    :rtype:


    """
    strucc,cvs_list = parse_pdb_and_generate_cvs(pdb_model, strucc=structure, peptide_length=peptide_length, one_model_per_nmr=True, only_reformat=only_reformat, process_only_chains=process_only_chains)
    return generate_matrix_and_graph(strucc, cvs_list, pdb_model, weight=weight, min_ah=min_ah, min_bs=min_bs, write_pdb=write_pdb,
                                     strictness_ah=strictness_ah, strictness_bs=strictness_bs, is_model=is_model, mixed_chains=True, path_base=path_base, verbose=verbose, process_only_chains=process_only_chains, interface=interface, domains=domains, compactness=compactness)

def parse_pdb_and_generate_cvs(pdbmodel, strucc=None, peptide_length=3, one_model_per_nmr=True, only_reformat=True, process_only_chains=None):
    if not strucc: strucc = Bioinformatics.get_structure("stru", pdbmodel)

    correct = []
    if not only_reformat:
        chains = Bio.PDB.Selection.unfold_entities(strucc, "C")
        seqs = [(chain.get_id(),"".join([Bioinformatics.AADICMAP[res.get_resname()] for res in Bio.PDB.Selection.unfold_entities(chain, "R") if res.has_id("CA") and res.has_id("C") and res.has_id("O") and res.has_id("N")])) for chain in chains]
        lich = []
        for s,seqr1 in enumerate(seqs):
            r1,seq1 = seqr1
            if r1 in lich: continue
            correct.append(r1)
            for t,seqr2 in enumerate(seqs):
                r2,seq2 = seqr2
                if t<=s or r2 in lich: continue
                z = galign(seq1,seq2)
                q = sum([1.0 for g in range(len(z[0][0])) if z[0][0][g] != "-" and z[0][1][g] != "-" and z[0][0][g]==z[0][1][g]])/len(z[0][0]) if len(z) > 0 else 0.0
                # print(z[0][0])
                # print(z[0][1])
                # print(z[0][2],q)
                # print()
                if q>=0.9: lich.append(r2)
        # for seq in seqs:
        #     print(seq)
        # print(correct)
    if process_only_chains:
        correct = list(process_only_chains)

    cvs_global, sep_chains = get_cvs(strucc, length_fragment=peptide_length, one_model_per_nmr=one_model_per_nmr, process_only_chains=correct)
    cvs_list = format_and_remove_redundance(cvs_global, sep_chains, only_reformat=only_reformat)
    return strucc,cvs_list

#@SystemUtility.timing
def generate_matrix_and_graph(strucc, cvs_list, pdb_model, weight="distance_avg", min_ah=4, min_bs=3, write_pdb=True,
                              strictness_ah=0.45, strictness_bs=0.20, is_model=False, mixed_chains=True,
                              maximum_distance=None, maximum_distance_bs=None, just_diagonal_plus_one=False, path_base="./",
                              process_only_chains=None, verbose=True, interface=False, domains=False, compactness=False):
    """
    NOTE CM: this documentation does not correspond with the parameters

    Annotates a protein pdb file with CVs and builds a graph in which each node is a secondary structure element
    or a coil fragment and the edge connecting two nodes is the metric distance between these fragments.

    :param pdb_model: The pdb file to annotate
    :type pdb_model: io.TextIOWrapper
    :param weight: The weight scheme for computing edge parameters. ["distance_avg","distance_min","distance_max"]
    :type weight: str
    :param min_ah: Minimum number of residues for an alpha helix to be accepted
    :type min_ah: int
    :param min_bs: Minimum number of residues for a beta strand to be accepted
    :type min_bs: int
    :param write_pdb: Write an annotated pdb file
    :type write_pdb: bool
    :param strictness_ah: strictness parameter threshold for accepting ah CVs
    :type strictness_ah: float
    :param strictness_bs: strictness parameter threshold for accepting bs CVs
    :type strictness_bs: float
    :param peptide_length: Define the peptide length for computing a CV
    :type peptide_length: int
    :return (g,strucc): The graph and the structure object representing the annotated pdb model
    :rtype (g,strucc): (igraph.Graph,Bio.PDB.Structure)
    """
    matrix, cvs_list, highd = get_3d_cvs_matrix(cvs_list, is_model, mixed_chains=mixed_chains,
                                                maximum_distance=maximum_distance,
                                                maximum_distance_bs=maximum_distance_bs,
                                                just_diagonal_plus_one=just_diagonal_plus_one)

    g = aleph_secstr(strucc, cvs_list, matrix, min_ah=min_ah, min_bs=min_bs, strictness_ah=strictness_ah, strictness_bs=strictness_bs, interface=interface)
    resil = [tuple(resi[1:4]) for frag in g.vs for resi in frag["reslist"]]
    toadd = []
    for r in Bioinformatics.get_list_of_residues(strucc):
        if process_only_chains is not None and r.get_full_id()[2] not in process_only_chains: continue
        if tuple(r.get_full_id()[1:4]) not in resil and r.get_resname() in Bioinformatics.AADICMAP:
            toadd.append({"reslist":[r.get_full_id()], "sstype":"coil", "resIdList":[r.get_full_id()], "cvids":[],"strictnesses":[], "cvls":[], "sequence":"".join([Bioinformatics.AADICMAP[r.get_resname()]])})
    s = igraph.Graph.Full(g.vcount()+len(toadd))
    i = 0
    for i, fr in enumerate(g.vs):
        for m in fr.attributes().keys():
            s.vs[i][m] = fr[m]
    i = g.vcount()
    for fr in toadd:
        s.vs[i]["reslist"] = fr["reslist"]
        s.vs[i]["sstype"] = fr["sstype"]
        s.vs[i]["resIdList"] = fr["resIdList"]
        s.vs[i]["cvids"] = fr["cvids"]
        s.vs[i]["cvls"] = fr["cvls"]
        s.vs[i]["sequence"] = fr["sequence"]
        s.vs[i]["strictnesses"] = fr["strictnesses"]
        #s.vs[i]["resforgraph"] = fr["resforgraph"]
        i += 1
    g = s

    g.vs["sheet"] = [None for _ in g.vs]

    g = aleph_internal_terstr(strucc, g, matrix, cvs_list)

    if domains: #To weight all edges equally
        for i,frag in enumerate(g.vs()):
            if frag['sstype']=='bs':
                g.vs[i]['sstype']='ah'

    g = aleph_terstr(strucc, g, matrix, cvs_list, weight=weight, verbose=verbose)

    if compactness: #To filter the graph and remove unconnected fragments
        g = g.vs.select(sstype_in=["ah", "bs"]).subgraph()
        nodes_to_delete = []

        for vertex in g.vs():
            unconnected = True
            for edge in vertex.all_edges():
                if edge['min']  < 20:
                    unconnected=False
            if unconnected==True:
                nodes_to_delete.append(vertex.index)

        g.delete_vertices(nodes_to_delete)

        #Alternative way of filtering but removes also distant edges and it can reaised a error in other functions
        # edges_filter = g.es.select(min_le=20)
        # g=g.subgraph_edges(edges_filter, delete_vertices=True)

    listar = [tuple(res[:-1]) for frag in g.vs for res in frag["reslist"] if frag["sstype"] in ["ah", "bs"]]

    #Calculating the percentage of secondary structure
    try:
        list_total_atoms = Bioinformatics.get_atoms_list(pdb_model)
        num = 0
        list_CA = [atom for atom in list_total_atoms if atom[2]=='CA']
        list_CA_no_disordered = []
        for atom in list_CA:
            if atom[5] == num:
                continue
            num = atom[5]
            list_CA_no_disordered.append(atom)
        n_res = len(list_CA_no_disordered)
        ah_content, bs_content = Bioinformatics.percentage_of_secondary_structure(g, n_res)
        dic_json = {'secondary_structure_content': {'ah':ah_content, 'bs':bs_content, 'number_total_residues': n_res}}
        modify_json_file('output.json', dic_json, annotation=True)
    except:
        print("Secondary structure content has not been calculated")

    pdbsearchin = Bioinformatics.generate_secondary_structure_record(g)
    #ss_record_dict = Bioinformatics3.write_csv_node_file(g, pdb_model,pdbsearchin, strictness_ah, strictness_bs)
    #Bioinformatics3.write_csv_edge_file(g, ss_record_dict)
    for model in strucc.get_list():
        reference = [atm for atm in Bio.PDB.Selection.unfold_entities(strucc[model.get_id()], "A") if (len(listar) == 0 or atm.get_parent().get_full_id() in listar)]
        pdbmod, cnv = Bioinformatics.get_pdb_from_list_of_atoms(reference, renumber=False, uniqueChain=False)
        pdbsearchin += pdbmod

    if write_pdb:
        fds = open(os.path.join(path_base, os.path.basename(pdb_model)[:-4] + "_input_search.pdb"), "w")
        fds.write(pdbsearchin)
        fds.close()
    #toadd = []

    return g, strucc, matrix, cvs_list, highd

def parse_hhr_file (pdb_model, g, write_pdb, hhr_file, path_base='./'):
    """Parsing hhr file and writing pdb containing alignment information in b-factors

    :param pdb_model: string containing the path to the target pdb model
    :type pdb_model: str
    :param g: graph containing the CV and their properties of the target model
    :type g: igraph.Graph
    :param write_pdb:
    :type: write_pdb: bool
    :param path_base: path where the output will be written
    :type path_base: str
    :return: g: graph containing the CV and their properties of the target model
    :rtype g: igraph.Graph
    """
    start_time_align = time.time()
    hhpred_to_bfac = {'|': 20, '+': 30, '.': 40, '-': 50, '=': 60, ' ': 60}
    hhpred_to_score = {'|': 1, '+': 0.75, '.': 0.5, '-': 0.25, '=': 0, ' ': 0}
    dict_sequence = Bioinformatics.get_sequence_from_pdb(pdb_model)
    strucc = Bioinformatics.get_structure("stru", pdb_model)

    if not os.path.isfile(hhr_file) or not hhr_file.endswith('.hhr'):
        print("ERROR: argument hhr_file should be a path to an existing .hhr file")
        raise Exception("argument hhr_file should be a path to an existing .hhr file")

    dict_homologs = Bioinformatics.reading_HHPRED_file(hhr_file)

    if not bool (dict_homologs):
        print("The .hhr file could not be read. File containing the alignment will not be written")
        return(g)

    dict_related = {}
    some_homolog = False

    for homolog in dict_homologs:
        for related in dict_homologs[homolog]['related']:
            dict_related[related] = homolog

    for model in strucc.get_list():
        for i, chain in enumerate(Bio.PDB.Selection.unfold_entities(strucc[model.get_id()], "C")):
            chain = str(chain).split('=')[1][:-1]
            modelo = (os.path.basename(pdb_model)[:-4] + '_' + chain).upper()
            if chain not in dict_sequence:
                continue
            start_aa = dict_sequence[chain]['first_res']
            align_sequence_1 = dict_sequence[chain]['sequence']


            for key in dict_homologs:
                if key.startswith(modelo[:-2]):
                    align_sequence_2 = ''.join(dict_homologs[key]['sequence'])
                    alignment = ''.join(dict_homologs[key]['alignment'])
                    some_homolog = True

            # if modelo in dict_homologs:
            #     align_sequence_2 = ''.join(dict_homologs[modelo]['sequence'])
            #     alignment = ''.join(dict_homologs[modelo]['alignment'])
            #     some_homolog = True
            # elif modelo in dict_related:
            #     related = dict_related[modelo]
            #     align_sequence_2 = ''.join(dict_homologs[related]['sequence'])
            #     alignment = ''.join(dict_homologs[related]['alignment'])
            #     some_homolog = True

            if some_homolog == False:
                print('The model %s has not been found in the alignment file' %(modelo))
                continue

            list_align = galign(align_sequence_1, align_sequence_2)

            for i, char in enumerate(list_align[0][0]):
                if char != '-':
                    index_1 = i
                    break

            for i, char in enumerate(list_align[0][1]):
                if char != '-':
                    index_2 = i
                    break

            conversion = - start_aa + index_1 - index_2

            for i, frag in enumerate(g.vs):
                frag['reslist_aligned'] = [resi for resi in frag['reslist'] if (
                        (resi[3][1] + conversion) >= 0 and (resi[3][1] + conversion) < len(alignment))]
                frag_alignment = [alignment[resi[3][1] + conversion] for resi in frag['reslist'] if (
                        (resi[3][1] + conversion) >= 0 and (resi[3][1] + conversion) < len(alignment))]

                frag['alignment_symbol'] = frag_alignment
                frag_alignment_bfac = [hhpred_to_bfac[symbol] for symbol in frag_alignment]
                frag['alignment_bfac'] = frag_alignment_bfac
                frag_alignment_score = [hhpred_to_score[symbol] for symbol in frag_alignment]
                frag['alignment_score'] = frag_alignment_score

        if some_homolog == False:
            continue

        reference = []
        for atm in Bio.PDB.Selection.unfold_entities(strucc[model.get_id()], "A"):
            set1 = set(atm.get_full_id()[0:-1])
            encontrado = False
            for frag in g.vs:
                if not 'reslist_aligned' in frag.attributes().keys() or frag['reslist_aligned'] == None:
                    continue
                for e, residue in enumerate(frag['reslist_aligned']):
                    set2 = set(residue)
                    if set1.issubset(set2):
                        try:
                            atm.set_bfactor(frag['alignment_bfac'][e])
                            reference.append(atm)
                            encontrado = True
                        except:
                            atm.set_bfactor(hhpred_to_bfac[' '])
                            reference.append(atm)
            if not encontrado:
                atm.set_bfactor(hhpred_to_bfac[' '])
                reference.append(atm)

        pdbmod_1, cnv_1 = Bioinformatics.get_pdb_from_list_of_atoms(reference, renumber=False,
                                                                    uniqueChain=False)
        pdbsearchin_1 = Bioinformatics.generate_secondary_structure_record(g)
        pdbsearchin_1 += pdbmod_1

        if write_pdb:
            fds = open(os.path.join(path_base, os.path.basename(pdb_model)[:-4] + "_alignment.pdb"), "w")
            fds.write(pdbsearchin_1)
            fds.close()

    if some_homolog == False:
        print("File containing the alignment will not be written")

    print("Time elapsed in alignment: {:.2f}s".format(time.time() - start_time_align))
    return(g)

def get_dictionary_from_community_clusters(graph, vclust, structure, writePDB=False, outputpath=None, header="",
                                           returnPDB=False, adding_t=0):
    global list_ids

    dict_res = {}
    pdbsearchin = ""
    if (writePDB and outputpath is not None) or returnPDB:
        pdbsearchin = header + "\n"

    for t, clust in enumerate(vclust):  # for every group in the community clustered graph
        listar = []
        for vertex in graph.vs:  # for every secondary structure element (vertex) in the clustered graph
            if vertex.index in clust:
                listar += [tuple(res[:-1]) for res in vertex["reslist"]]

        for model in structure.get_list():
            reference = []
            for chain in model.get_list():
                for residue in chain.get_list():
                    if len(listar) == 0 or residue.get_full_id() in listar:
                        reference += residue.get_unpacked_list()
                        dict_res[residue.get_full_id()] = 'group' + str(t + adding_t)
            if (writePDB and outputpath is not None) or returnPDB:
                pdbmod, cnv = Bioinformatics.get_pdb_from_list_of_atoms(reference, renumber=True, uniqueChain=True,
                                                                        chainId=list_ids[t + adding_t])
                pdbsearchin += pdbmod

        if writePDB and outputpath is not None:
            fds = open(outputpath, "w")
            fds.write(pdbsearchin)
            fds.close()

    if not returnPDB:
        return dict_res
    else:
        return dict_res, pdbsearchin

def get_community_clusters(graph, algorithm="fastgreedy", n=None, print_dendo=None, return_dendo=False, weight="weight"):
    vdendo = None
    if algorithm.lower() == "fastgreedy":
        vdendo = graph.community_fastgreedy(weights=weight)
        try:
            vclust = vdendo.as_clustering(n=n)
        except:
            vclust = None

        if print_dendo is not None:
            try:
                igraph.drawing.plot(vdendo, target=print_dendo, bbox=(0, 0, 800, 800))
            except:
                try:
                    igraph.drawing.plot(vdendo, target=print_dendo[:-4] + ".svg", bbox=(0, 0, 800, 800))
                except:
                    pass
    elif algorithm.lower() == "infomap":
        vclust = graph.community_infomap(edge_weights=weight)
    elif algorithm.lower() == "eigenvectors":
        vclust = graph.community_leading_eigenvector(weights=weight)
    elif algorithm.lower() == "label_propagation":
        vclust = graph.community_label_propagation(weights=weight)
    elif algorithm.lower() == "community_multilevel":
        vclust = graph.community_multilevel(weights=weight, return_levels=False)
    elif algorithm.lower() == "edge_betweenness":
        vdendo = graph.community_edge_betweenness(directed=False, weights=weight)
        try:
            vclust = vdendo.as_clustering(n=n)
        except:
            vclust = None

        if print_dendo is not None:
            try:
                igraph.drawing.plot(vdendo, target=print_dendo, bbox=(0, 0, 800, 800))
            except:
                try:
                    igraph.drawing.plot(vdendo, target=print_dendo[:-4] + ".svg", bbox=(0, 0, 800, 800))
                except:
                    pass
    elif algorithm.lower() == "spinglass":
        vclust = graph.community_spinglass(weights=weight)
    elif algorithm.lower() == "walktrap":
        vdendo = graph.community_walktrap(weights=weight)
        if print_dendo is not None:
            try:
                igraph.drawing.plot(vdendo, target=print_dendo, bbox=(0, 0, 800, 800))
            except:
                try:
                    igraph.drawing.plot(vdendo, target=print_dendo[:-4] + ".svg", bbox=(0, 0, 800, 800))
                except:
                    pass
        try:
            vclust = vdendo.as_clustering(n=n)
        except:
            vclust = None
    else:
        vclust = None
    if not return_dendo:
        return vclust
    else:
        return vclust, vdendo


def pack_beta(vclust,graph,pack_beta_sheet):
    if not pack_beta_sheet:
        return vclust

    y = vclust.subgraphs()
    while 1:
        merge = False
        indexes = []
        r = None
        for i, clu1 in enumerate(y):
            for j, clu2 in enumerate(y):
                if j > i:
                    # print([t["sstype"] for t in clu1.vs],[t["sstype"] for t in clu2.vs])
                    # print([t["sstype"]=="bs" for t in clu1.vs],[t["sstype"]=="bs" for t in clu2.vs],[s["mean"] for s in sorted(graph.es.select(lambda x: (graph.vs[x.source]["name"] in clu1.vs["name"] and graph.vs[x.target]["name"] in clu2.vs["name"]) or (graph.vs[x.source]["name"] in clu2.vs["name"] and graph.vs[x.target]["name"] in clu1.vs["name"])), key=lambda e: e["avg"])])
                    # print([(t["sequence"],t["reslist"][0][2]) for t in clu1.vs],[(t["sequence"],t["reslist"][0][2]) for t in clu2.vs])
                    if all([t["sstype"] == "bs" for t in clu1.vs]) and all([t["sstype"] == "bs" for t in clu2.vs]):
                        if len(set(clu1.vs["sheet"])&set(clu2.vs["sheet"]))==1:
                            r = graph.vs.select(name_in=clu1.vs["name"] + clu2.vs["name"]).subgraph()
                            merge = True
                            indexes.append(i)
                            indexes.append(j)
                            break
            if merge:
                break
        if merge:
            # print("merging",indexes)
            z = [q for p, q in enumerate(y) if p not in indexes] + [r]
            y = z
        else:
            break

    # ul = []
    # for o in graph.vs["name"]:
    #     for l, p in enumerate(y):
    #         if o["name"] in p.vs["name"]:
    #             ul.append(l)
    #             break

    ul = [l for o in graph.vs["name"] for l, p in enumerate(y) if o in p.vs["name"]]

    # print(graph.vcount())
    # print(ul)
    vclust = igraph.VertexClustering(graph, membership=ul)
    return vclust


def get_community_clusters_one_step(algo, graph_input_d, structure, pdb_search_in, pathpdb, n=None, print_dendo=None,
                                    return_dendo=False, use_dendo=None, use_spantree=True, write_pdb=True,
                                    weight="weight", pack_beta_sheet=False, homogeneity=False, write_graph_img = False, path_base="./"):


    # for u in sorted(graph_input.es["weight"],reverse=True):
    #     print("before",u)
    # print()
    graph_input = graph_input_d.copy()
    graph_input = graph_input.vs.select(sstype_in=["ah", "bs"]).subgraph()
    graph_input.vs["name"] = [str(vertex["reslist"][0][3][1]) + "_" + str(vertex["reslist"][0][2]) for vertex in graph_input.vs]


    if pack_beta_sheet:
        graph_input.es["weight"] = [100*(100.0/e["avg"]) if graph_input.vs[e.source]["sstype"] == graph_input.vs[e.target]["sstype"] else e["weight"] for e in graph_input.es]
        graph_input.es["spantree_weight"] = [e["avg"]/100.0 if graph_input.vs[e.source]["sstype"] == graph_input.vs[e.target]["sstype"] else e["spantree_weight"] for e in graph_input.es]


    if use_spantree:
        graph = graph_input.spanning_tree(weights="spantree_weight")
    else:
        graph = graph_input
    # graph.write_graphml(os.path.join("./", os.path.basename(pathpdb)[:-4] + "_spantree.graphml"))
    # for u in sorted(graph.es["weight"],reverse=True):
    #     print("after",u)
    # print()

    if not homogeneity:
        if not return_dendo and use_dendo is None:
            vclust = get_community_clusters(graph, algorithm=algo, n=n, print_dendo=print_dendo, return_dendo=return_dendo,  weight=weight)
        elif not return_dendo:
            vclust = use_dendo.as_clustering(n=n)
        else:
            vclust, vdendo = get_community_clusters(graph, algorithm=algo, n=n, print_dendo=print_dendo,
                                                    return_dendo=return_dendo, weight=weight)
    else:
        return_dendo = True
        vclust, vdendo = get_community_clusters(graph, algorithm=algo, n=n, print_dendo=print_dendo, return_dendo=return_dendo, weight=weight)
        l = [pack_beta(vdendo.as_clustering(n=i),graph,pack_beta_sheet) for i in range(1,graph.vcount())]
        best_i = None
        best_score = -1.0
        for w,cluster in enumerate(l):
            q=tuple([tuple([t["sstype"] for t in c.vs]) for c in cluster.subgraphs()])
            if not pack_beta_sheet:
                sample = [c.vcount() for c in cluster.subgraphs()]
            else:
                sample = [c.vcount() for c in cluster.subgraphs() if len(set(tuple(c.vs["sstype"])))==1 and c.vs["sstype"][0]=="ah"]

            if len(sample) == 1 or (len(sample) == 0 and pack_beta_sheet):
                uq = 0.0
                CV = 0.0
                normCV = 0.1
            elif len(sample) == 0:
                print("Error!!!!, Clustering has not grouped anything?")
                raise Exception("Error!!!!, Clustering has not grouped anything?")
            else:
                uq = numpy.mean(sample)
                sq = numpy.std(sample)
                CV = sq / uq
                normCV = (((CV - 0) * (1 - 0)) / (numpy.sqrt(len(sample) - 1) - 0)) + 0

            n = graph.vcount()
            # (((OldValue - OldMin) * (NewMax - NewMin)) / (OldMax - OldMin)) + NewMin

            corrected_score = cluster.modularity+(uq/n)-(normCV) #cluster.modularity*uq*(len(q)-len(set(q)))
            if  corrected_score > best_score:
                best_i = w
                best_score = corrected_score
            #print(w,cluster.modularity,q,len(q))
            #print(w,cluster.modularity,set(q),len(set(q)))
            #print(w,cluster.modularity,uq/n,normCV,corrected_score)

        if best_i is None:
            return None

        print("Best cut",best_i+1,"modularity",best_score)
        vclust = l[best_i]

    if pack_beta_sheet:
        vclust = pack_beta(vclust,graph,pack_beta_sheet)

    if write_pdb:
        get_dictionary_from_community_clusters(graph, vclust, structure, writePDB=True, outputpath=os.path.join(path_base,
                                                                                                                os.path.basename(
                                                                                                                    pathpdb)[
                                                                                                                :-4] + "_" + algo + "_distclust.pdb"),
                                               header=pdb_search_in)
    layout = graph.layout("kk")

    if write_graph_img:
        try:
            write_image_from_graph(vclust,os.path.join(path_base, os.path.basename(pathpdb)[:-4] + "_" + algo + "_clustering.png"))
        except:
            # vclust.write_graphml(os.path.join("./", os.path.basename(pathpdb)[:-4] + "_clustering.graphml"))
            pass


    dic_json = {pathpdb: {'number of groups': len(vclust), 'size of groups': list(len(elem) for elem in vclust)}}
    modify_json_file('output.json', dic_json, decomposition=True)

    if return_dendo:
        return vclust, vdendo
    else:
        return vclust

@SystemUtility.timing
def generate_graph_for_cvs(graph_full, matrix, cvs_list, peptide_length=3):
    # tuplone = Bioinformatics.getFragmentListFromPDBUsingAllAtoms(pdb_model, False)
    # strucc = tuplone[0]

    if isinstance(cvs_list[0][0], list):
        cvs_list = [lis for l in cvs_list for lis in l]

    sstype_map = {"ah": 1, "bs": 2, "coil": 3}

    g = igraph.Graph()
    g.add_vertices(len(cvs_list))

    for z, cvs in enumerate(cvs_list):
        g.vs[z]["secstr"] = None
        for p, frag in enumerate(graph_full.vs):
            a = tuple(tuple(x) for x in cvs[4])
            b = tuple(tuple(x) for x in frag["reslist"])
            if set(a).issubset(b):
                g.vs[z]["secstr"] = p
                ul = get_vseq_neighbours_fragments(graph_full, frag, sortmode="avg")
                if len(ul) > 0:
                    g.vs[z]["nearsecstr"] = ul[0].index
                else:
                    g.vs[z]["nearsecstr"] = None
                #g.vs[z]["com"] = frag["com"]
                g.vs[z]["unique_cv"] = frag["unique_cv"]
                g.vs[z]["cvs"] = z
                g.vs[z]["-1"] = z - 1 if (z - 1 > 0 and g.vs[z - 1]["secstr"] == p) else None
                g.vs[z]["+1"] = None
                if (z - 1 > 0 and g.vs[z - 1]["secstr"] == p):
                    g.vs[z - 1]["+1"] = z

                g.vs[z]["reslist"] = cvs[4]
                g.vs[z]["sstype"] = frag["sstype"]
                g.vs[z]["sequence"] = "".join(map(lambda x: Bioinformatics.AADICMAP[x[4]], cvs[4]))
                g.vs[z]["sequence_secstr"] = frag["sequence"]
                g.vs[z]["isStart"] = True if g.vs[z]["reslist"] == frag["reslist"][:peptide_length] else False
                g.vs[z]["isMiddle"] = True if g.vs[z]["reslist"] == frag["reslist"][int(len(frag["reslist"]) / 2):int(
                    len(frag["reslist"]) / 2) + peptide_length] else False
                g.vs[z]["isEnd"] = True if g.vs[z]["reslist"] == frag["reslist"][
                                                                 len(frag["reslist"]) - peptide_length:] else False
                g.vs[z]["isSpecial"] = True if g.vs[z]["isStart"] or g.vs[z]["isMiddle"] or g.vs[z]["isEnd"] else False
                g.vs[z]["posandsec"] = ((frag["reslist"].index(g.vs[z]["reslist"][0]) + 1) * 100) + p
                g.vs[z]["pos"] = (frag["reslist"].index(g.vs[z]["reslist"][0]) + 1)
                break

    g = g.vs.select(lambda vertex: vertex["secstr"] is not None).subgraph()

    for frag1 in g.vs:
        for frag2 in g.vs.select(lambda vertex: vertex.index > frag1.index):
            g.add_edge(frag1.index, frag2.index)
            value1 = matrix[(frag1["cvs"], frag2["cvs"])]
            value2 = compute_instruction(frag1["unique_cv"], frag2["unique_cv"], unique_fragment_cv=True)
            z = [value1[2], value1[3], value1[4], sstype_map[frag1["sstype"]], sstype_map[frag2["sstype"]]]
            r = [value2[2], value2[3], value2[4]]
            w = [compute_instruction(cvs_list[frag1["cvs"]],frag2["unique_cv"])[3], compute_instruction(cvs_list[frag2["cvs"]],frag1["unique_cv"])[3]]
            g.es[g.get_eid(frag1.index, frag2.index)]["weight"] = 100.0 / value1[3]
            g.es[g.get_eid(frag1.index, frag2.index)]["metric"] = z+r+w+[len(set([frag1["secstr"], frag2["secstr"]]))]#+[frag1["pos"], frag2["pos"]]
            g.es[g.get_eid(frag1.index, frag2.index)]["distance"] = value1[3]

    # for frag1 in g.vs:
    #     for frag2 in g.vs.select(lambda vertex: vertex.index > frag1.index):
    #         g.es[g.get_eid(frag1.index, frag2.index)]["metric"] += [g.es.select(lambda x: x["distance"] >=5.0 <)]
    g.vs["name"] = [str(vertex.index) + "_" + str(vertex["reslist"][0][2]) for vertex in g.vs]

    return g


@SystemUtility.timing
def decompose_by_community_clustering(reference, strictness_ah, strictness_bs, peptide_length, pack_beta_sheet,
                                      max_ah_dist, min_ah_dist, max_bs_dist, min_bs_dist, write_graphml, write_pdb,
                                      homogeneity, algorithm="fastgreedy", hhr_file='', write_graph_img=False):
    # NOTE: PARAMETRIZATION=======
    write_graphml = bool(write_graphml)
    sort_mode = "avg"

    if strictness_ah < 0:
        strictness_ah = 0.45
    if strictness_bs < 0:
        strictness_bs = 0.20
    if min_bs_dist < 0:
        min_bs_dist = 0
    if max_bs_dist < 0:
        max_bs_dist = 15
    if min_ah_dist < 0:
        min_ah_dist = 0
    if max_ah_dist < 0:
        max_ah_dist = 20

    weight = "distance_avg"
    pep_len = int(peptide_length)

    pdb_search_in = ""
    f = open(reference, "r")
    all_lines = f.readlines()
    f.close()
    for line in all_lines:
        if line[:6] in ["TITLE ", "CRYST1", "SCALE1", "SCALE2", "SCALE3"]:
            pdb_search_in += line

    graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
        reference, weight=weight, strictness_ah=strictness_ah, strictness_bs=strictness_bs, peptide_length=pep_len,
        write_pdb=write_pdb)

    if hhr_file:
        graph_full_reference = parse_hhr_file(reference, graph_full_reference, write_pdb, hhr_file)

    graph_full_reference.vs["name"] = [str(vertex["reslist"][0][3][1]) + "_" + str(vertex["reslist"][0][2]) for vertex
                                       in graph_full_reference.vs]
    graph_no_coil_reference = graph_full_reference.vs.select(sstype_in=["ah", "bs"]).subgraph()
    eigen_ref = graph_no_coil_reference.evcent(directed=False, scale=True, weights=graph_no_coil_reference.es["weight"],
                                               return_eigenvalue=True)


    # print("EIGEN Centrality values:", eigen_ref, len(eigen_ref[0]), len(graph_no_coil_reference.vs))
    graph_no_coil_reference.vs["eigen"] = eigen_ref[0]
    graph_no_coil_reference.vs["Label"] = [d["sstype"] + "_" + d["name"] for t, d in
                                           enumerate(graph_no_coil_reference.vs)]

    if write_graphml:
        graph_no_coil_reference.write_graphml(os.path.join("./", os.path.basename(reference)[:-4] + "_graphref.graphml"))

    for i, frag in enumerate(graph_no_coil_reference.vs):
        print(
            frag.index, "**", frag["sstype"], frag["reslist"][0][2], frag["reslist"][0][3][1], "--",
            frag["reslist"][-1][3][1], frag["sequence"], "CVIDS: ", frag["cvids"][0], frag["cvids"][
                int(len(frag["cvids"]) / 2)], frag["cvids"][-1], "SHEET ID: "+str(frag["sheet"]) if frag["sstype"] == "bs" else "")
        if sort_mode == "avg":
            print("\t\tSORTED LIST OF FRAGMENT BY AVG DISTANCE:")
            print("\t\tPOS:\tN_FRAG:\tAVG_DIST:")
        else:
            print("\t\tSORTED LIST OF FRAGMENT BY MIN DISTANCE:")
            print("\t\tPOS:\tN_FRAG:\tMIN_DIST:")

        for j, edge2 in enumerate(get_eseq_neighbours_fragments(graph_no_coil_reference, frag, sortmode=sort_mode)):
            print("\t\t" + str(j) + "\t" + str(get_connected_fragment_to_edge(frag, edge2).index) + "\t" + str(
                edge2[sort_mode]))

    print("EXECUTING COMMUNITY CLUSTERING....")

    get_community_clusters_one_step(algorithm, graph_no_coil_reference, structure_reference, pdb_search_in, os.path.basename(reference),pack_beta_sheet=pack_beta_sheet,homogeneity=homogeneity, write_graph_img=write_graph_img)


@SystemUtility.timing
def annotate_pdb_model(reference, strictness_ah, strictness_bs, peptide_length, width_pic, height_pic, write_graphml,
                       write_pdb, hhr_file = None, verbose=False, compactness=False):


    # NOTE: PARAMETRIZATION=======
    # write_graphml = bool(write_graphml)
    sort_mode = "avg"

    if strictness_ah < 0:
        strictness_ah = 0.45
    if strictness_bs < 0:
        strictness_bs = 0.20

    weight = "distance_avg"
    pep_len = int(peptide_length)

    pdb_search_in = ""
    f = open(reference, "r")
    all_lines = f.readlines()
    f.close()
    for line in all_lines:
        if line[:6] in ["TITLE ", "CRYST1", "SCALE1", "SCALE2", "SCALE3"]:
            pdb_search_in += line

    graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
        reference, weight=weight, strictness_ah=strictness_ah, strictness_bs=strictness_bs, peptide_length=pep_len,
        write_pdb=write_pdb, compactness=compactness)

    if hhr_file:
        graph_full_reference = parse_hhr_file(reference, graph_full_reference, write_pdb, hhr_file)

    graph_full_reference.vs["name"] = [str(vertex["reslist"][0][3][1]) + "_" + str(vertex["reslist"][0][2]) for vertex in graph_full_reference.vs]

    if write_graphml:
        graph_full_reference.write_graphml(os.path.join("./", os.path.basename(reference)[:-4] + "_graphref.graphml"))
        a = list(range(len(cvs_reference)))
        q = []
        for frag in graph_full_reference.vs:
            for b in a:
                if int(b) in frag["cvids"] or b in frag["cvids"]:
                    q.append((frag["sstype"],frag["name"]))

        newg = igraph.Graph(n=len(a))
        newg.vs["name"] = list(a)
        newg.vs["sstype"] = [z[0] for z in q]
        newg.vs["fragid"] = [z[1] for z in q]

        s = [(newg.vs.find(name=key1).index, newg.vs.find(name=key2).index, matrix_reference[tuple(sorted([cvs_reference[key1][0], cvs_reference[key2][0]]))][2], matrix_reference[tuple(sorted([cvs_reference[key1][0], cvs_reference[key2][0]]))][3]) for key1,ele1 in enumerate(cvs_reference) for key2,ele2, in enumerate(cvs_reference)]
        f = [(b[0],b[1]) for b in s]
        j = [b[2] for b in s]
        g = [b[3] for b in s]
        newg.add_edges(f)
        newg.es["distance"] = g
        newg.es["angle"] = j

        newg.write_graphml(os.path.join("./", os.path.basename(reference)[:-4] + "_cvs.graphml"))

    #NOTE: Writing here a pdb annotating the strictnesses in bfactors:
    qd = {}
    EXTREME_COIL = 0.0
    for i, frag in enumerate(graph_full_reference.vs):
        # print(i)
        # print(frag["sstype"])
        # print(frag["cvids"], len(frag["cvids"]))
        # print(frag["strictnesses"], len(frag["strictnesses"]))
        # print(frag["sequence"],len(frag["sequence"]))
        # print(frag["reslist"], len(frag["reslist"]))
        assert(len(frag["cvids"]) == 0 or len(frag["cvids"]) == len(frag["strictnesses"]))
        if len(frag["reslist"]) == len(frag["strictnesses"]):
            qd.update({tuple(t[:-1]):frag["strictnesses"][p] for p,t in enumerate(frag["reslist"])})
        elif len(frag["cvids"]) == 0 and len(frag["reslist"]) > 0:
            qd.update({tuple(t[:-1]):EXTREME_COIL for t in frag["reslist"]})
        else:
            qd.update({tuple(t[:-1]):frag["strictnesses"][int(p/pep_len)] for p,t in enumerate(frag["reslist"])})

    qc = {}
    for res in Bioinformatics.get_list_of_residues(structure_reference):
        for atm in res:
            qc[tuple(atm.get_full_id())] = float(atm.get_bfactor())
            # print(res.get_full_id())
            # print(sorted(qd.keys()))
            # print()
            if tuple(res.get_full_id()) in qd:
                atm.set_bfactor(qd[tuple(res.get_full_id())])
            else:
                atm.set_bfactor(EXTREME_COIL)

    stren = Bioinformatics.get_pdb_from_structure(structure_reference, None)

    with open("strictnesses.pdb","w") as f: f.write(stren)

    for res in Bioinformatics.get_list_of_residues(structure_reference):
        for atm in res:
            atm.set_bfactor(qc[tuple(atm.get_full_id())])


    graph_no_coil_reference = graph_full_reference.vs.select(sstype_in=["ah", "bs"]).subgraph()


    try:
        eigen_ref = graph_no_coil_reference.evcent(directed=False, scale=True, weights=graph_no_coil_reference.es["weight"],
                                               return_eigenvalue=True)
        # print("EIGEN Centrality values:", eigen_ref, len(eigen_ref[0]), len(graph_no_coil_reference.vs))
        graph_no_coil_reference.vs["eigen"] = eigen_ref[0]
    except:
        graph_no_coil_reference.vs["eigen"] = [0 for x in graph_no_coil_reference.vs]

    graph_no_coil_reference.vs["Label"] = [d["sstype"] + "_" + d["name"] for t, d in enumerate(graph_no_coil_reference.vs)]
    # graph_no_coil_reference.write_graphml(
    #     os.path.join("./", os.path.basename(reference)[:-4] + "_graphref.graphml"))

    pdbsearchin = ""


    for i, frag in enumerate(graph_no_coil_reference.vs):
        print(
            frag.index, "**", frag["sstype"], frag["reslist"][0][2], frag["reslist"][0][3][1], "--",
            frag["reslist"][-1][3][1], frag["sequence"], "CVIDS: ", frag["cvids"][0], frag["cvids"][
                int(len(frag["cvids"]) / 2)], frag["cvids"][-1])
        if sort_mode == "avg":
            print("\t\tSORTED LIST OF FRAGMENT BY AVG DISTANCE:")
            print("\t\tPOS:\tN_FRAG:\tAVG_DIST:")
        else:
            print("\t\tSORTED LIST OF FRAGMENT BY MIN DISTANCE:")
            print("\t\tPOS:\tN_FRAG:\tMIN_DIST:")

        for j, edge2 in enumerate(get_eseq_neighbours_fragments(graph_no_coil_reference, frag, sortmode=sort_mode)):
            print("\t\t" + str(j) + "\t" + str(get_connected_fragment_to_edge(frag, edge2).index) + "\t" + str(
                edge2[sort_mode]))

    if verbose:
        list_distance = [edge['avg'] for edge in graph_no_coil_reference.es()] #In case we want the max distance from the complete structure
        list_distance = sorted(list_distance, reverse=True)
        print("\nThe 2 largest distances between any kind of secondary structure fragment are", round(list_distance[0],1), round(list_distance[1],1))

        graph_no_coil_reference_ah = graph_full_reference.vs.select(sstype_in=["ah"]).subgraph() #Only for ah
        list_distance = [edge['avg'] for edge in  graph_no_coil_reference_ah.es()]
        list_distance = sorted(list_distance, reverse=True)
        print("The 2 largest distances between alpha helix fragments are", round(list_distance[0], 1), round(list_distance[1], 1))

        graph_no_coil_reference_bs = graph_full_reference.vs.select(sstype_in=["bs"]).subgraph() #Only for bs
        list_distance = [edge['avg'] for edge in  graph_no_coil_reference_bs.es()]
        list_distance = sorted(list_distance, reverse=True)
        print("The 2 largest distances between beta strand fragments are", round(list_distance[0], 1), round(list_distance[1], 1) , '\n')


    f = open("./CVs.txt", "w")
    dic_json = {}
    dic_cvstxt = {}

    for frag in graph_full_reference.vs:
        for element in frag['cvids']:
            dic_cvstxt[element]= {'sstype': frag['sstype'], 'start':'', 'end':''}

    for i, element in enumerate(cvs_reference):
        try:
            dic_cvstxt[element[0]]['start']=numpy.array(element[2]).tolist()
            dic_cvstxt[element[0]]['end'] = numpy.array(element[3]).tolist()
            f.write(dic_cvstxt[element[0]]['sstype'] + '\n' + str(numpy.array(element[2]).tolist()) + '\n' + str(numpy.array(element[3]).tolist()) + '\n')
        except:
            f.write('' + '\n' + str(numpy.array(element[2]).tolist()) + '\n' + str(numpy.array(element[3]).tolist()) + '\n')
        dic_json[i] = {'start': numpy.array(element[2]).tolist(), 'end': numpy.array(element[3]).tolist()}

    modify_json_file('output.json', dic_json, cvs=True)
    f.close()

    # NOTE: Graphs of the secondary structures
    plt.figure(figsize=(width_pic, height_pic))
    ax = plt.subplot()
    x = []
    y = []
    z = []
    q = []
    dizioq = {}
    cmap = {"ah": "r", "bs": "y", "coil": "b"}
    for frag in graph_full_reference.vs:
        x += frag["cvids"]
        y += frag["cvls"]
        z += [frag["sstype"]] * len(frag["cvls"])
        q += [str(ty[3][1]) + "_" + str(ty[2]) for ty in frag["reslist"]]
        # print(frag["cvids"],frag["cvls"],[frag["sstype"]]*len(frag["cvls"]),[str(ty[3][1])+"_"+str(ty[2]) for ty in frag["reslist"]])
        # print(len(frag["cvids"]),len(frag["cvls"]),len([frag["sstype"]]*len(frag["cvls"])),len([str(ty[3][1])+"_"+str(ty[2]) for ty in frag["reslist"][:-2]]))

    dic_json = {'CVL_pic': {'cvids': x, 'cvls':y, 'sstype': z, 'label': q}}
    modify_json_file('output.json', dic_json, plot=True)

    ax.scatter(x, y, s=75, color=[cmap[u] for u in z])
    ax.plot(x, y, color='green', linewidth=2)
    for xyq in zip(x, y, q):
        ax.annotate(xyq[2].replace("_", ""), xy=xyq[:2], textcoords='data', ha='center', va="bottom")
        dizioq[xyq[0]] = xyq[2]
    ax.axhline(y=2.2, color='b')
    ax.axhline(y=1.4, color='b')
    trans = transforms.blended_transform_factory(ax.transAxes, ax.transData)
    ax.annotate("alpha", xy=(1.01, 2.2), xycoords=trans, clip_on=False, va='center')
    ax.annotate("beta", xy=(1.01, 1.4), xycoords=trans, clip_on=False, va='center')
    plt.title('' + os.path.basename(reference) + ' CVLs', fontsize=20)
    plt.ylabel('CVL', fontsize=16)
    plt.xlabel('CV id', fontsize=16)
    # plt.xlim(numpy.arange(numpy.min(x), numpy.max(x), 1.0))
    plt.xticks(numpy.arange(numpy.min(x), numpy.max(x), 1.0))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Alpha', markerfacecolor='r', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Beta', markerfacecolor='y', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Coil', markerfacecolor='b', markersize=10)]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.savefig(os.path.join("./", os.path.basename(reference)[:-4] + "_CVLs.png"))
    plt.close()

    f = open(os.path.join("./", os.path.basename(reference)[:-4] + "_CVLs.txt"), "w")
    for opo, l in enumerate(x):
        f.write(str(x[opo]) + "\t" + str(y[opo]) + "\t" + str(q[opo]) + "\t" + str(z[opo]) + "\n")
    # for opo, l in enumerate(x):
    #     f.write(str(x[opo]) + "\t" + str(y[opo]) + "\t" + str(z[opo]) +"\n")
    f.close()


    # NOTE: Graph of the distances of the secondary structures
    plt.figure(figsize=(width_pic, height_pic))
    x = []
    y = []
    z = []
    q = []
    r = []
    p = []
    s = []
    cmap = {"ah": "r", "bs": "y", "coil": "b"}

    for frag in graph_full_reference.vs:
        if frag["sstype"] == "bs":
            for cvid in frag["cvids"]:
                for frag2 in graph_full_reference.vs:
                    if frag2.index != frag.index:
                        suball = [valo for valo in
                                  graph_full_reference.es[graph_full_reference.get_eid(frag.index, frag2.index)]["alls"]
                                  if cvid in valo[:2] and valo[3] <= 10.0]
                        x += [cvid] * len(suball)
                        y += [ty[3] for ty in suball]
                        s += [ty[2] for ty in suball]
                        q += [ty[0] if ty[0] != cvid else ty[1] for ty in suball]
                        z += [frag2["sstype"]] * len(suball)

    f = open(os.path.join("./", os.path.basename(reference)[:-4] + "_distances.txt"), "w")
    for opo, l in enumerate(x):
        f.write(str(x[opo]) + "\t" + str(q[opo]) + "\t" + str(y[opo]) + "\t" + str(s[opo]) + "\t" + str(
            z[opo]) + "\t" + str(dizioq[x[opo]]) + "\t" + str(dizioq[q[opo]]) + "\n")
        # f.write(str(x[opo])+"\t"+str(q[opo])+"\t"+str(y[opo])+"\t"+str(z[opo])+"\t"+"\n")
    f.close()

    # NOTE: Graph of the angles in secondary structures
    plt.figure(figsize=(width_pic, height_pic))
    x = []
    y = []
    z = [] #Secodnary structure
    q = [] #Cvs label
    o = [] #Fragment first CV id
    n = [] #Fragment angle mean
    l = [] #Fragmetn label
    k = [] #Ca distance
    cmap = {"ah": "r", "bs": "y", "coil": "b"}



    graph_full_reference = graph_full_reference.vs.select(lambda x: len(x["cvids"]) > 0).subgraph()

    for frag in graph_full_reference.vs:
        # print frag["cvids"]
        # print [ty[0] for ty in frag["alls"]]
        # print [frag["sstype"]] * len(frag["cvls"][:-1])
        # print
        # print
        if frag.index < len(graph_full_reference.vs) - 1:
            o += [frag["cvids"][0]]
            n += [graph_full_reference.es[graph_full_reference.get_eid(frag.index, frag.index + 1)]["mean"][0]]
            l += [str(frag["reslist"][0][3][1]) + "-" + str(frag["reslist"][-1][3][1]) + "_" + str(
                frag["reslist"][0][2]) + "--" + str(
                graph_full_reference.vs[frag.index + 1]["reslist"][0][3][1]) + "-" + str(
                graph_full_reference.vs[frag.index + 1]["reslist"][-1][3][1]) + "_" + str(
                graph_full_reference.vs[frag.index + 1]["reslist"][0][2])]

        if len(frag["alls"])!= len(frag["reslist"][:-2])-1:
            y += [ty[0] for ty in frag["alls"][2:]]
            k += [ty[1] for ty in frag["alls"][2:]]
            x += frag["cvids"][:-3]
            z += [frag["sstype"]] * len(frag["cvls"][:-3])
        else:
            x += frag["cvids"][:-1]
            y += [ty[0] for ty in frag["alls"]]
            z += [frag["sstype"]] * len(frag["cvls"][:-1])
            k += [ty[1] for ty in frag["alls"]]

        q += [str(frag["reslist"][ty][3][1]) + "_" + str(frag["reslist"][ty][2]) + "-" + str(
            frag["reslist"][ty + 1][3][1]) + "_" + str(frag["reslist"][ty + 1][2]) for ty in
              range(len(frag["reslist"][:-2]) - 1)]

        #print(str(frag["reslist"][3][1]))

    # print("LEN",len(x),len(y))
    # print(x)
    # print(y)

    dic_json = {'Angles_pic': {'cvids': x, 'cv_anles':y, 'sstype': z, 'label': q, 'ca_distance' : k, 'fragment_start' : o, 'fragment_mean' : n, 'fragment_label':l}}
    modify_json_file('output.json', dic_json, plot=True)

    plt.scatter(x, y, s=75, color=[cmap[u] for u in z])
    plt.plot(x, y, linewidth=2, color="g")
    plt.scatter(o, n, linewidth=2, color="g")

    f = open(os.path.join("./", os.path.basename(reference)[:-4] + "_Angles.txt"), "w")

    for xyq in zip(x, y, q, z):
        plt.annotate(xyq[2].replace("_",""), xy=xyq[:2], textcoords='data', ha='center', va="bottom")
        f.write(str(xyq[0]) + "\t" + str(xyq[1]) + "\t" + str(xyq[2]) + "\t" + str(xyq[3]) + "\n")
    # for xy in zip(x, y, z):
    #    f.write(str(xy[0]) + "\t" + str(xy[1]) + "\t" + str(xy[2]) +"\n")

    for onl in zip(o, n, l, z):
        plt.annotate(onl[2].replace("_",""), xy=onl[:2], textcoords='data', ha='center', va="bottom")
        # f.write(str(onl[0]) + "\t" + str(onl[1]) + "\t" + str(onl[2]) + "\t" + str(onl[3]) +"\n")

    # plt.axhline(y=2.2, color='b')
    # plt.axhline(y=1.4, color='b')
    plt.title('' + os.path.basename(reference) + ' Angles', fontsize=20)
    plt.ylabel('Angle', fontsize=16)
    plt.xlabel('CV id', fontsize=16)
    # # plt.xlim(numpy.arange(numpy.min(x), numpy.max(x), 1.0))
    plt.xticks(numpy.arange(numpy.min(x), numpy.max(x), 1.0))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Alpha', markerfacecolor='r', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Beta', markerfacecolor='y', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Coil', markerfacecolor='b', markersize=10)]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.savefig(os.path.join("./", os.path.basename(reference)[:-4] + "_Angles.png"))
    plt.close()

    f.close()

    #NOTE: Graph of the Ca-Ca distances
    plt.figure(figsize=(width_pic, height_pic))
    plt.scatter(x, k, s=75, color=[cmap[u] for u in z])
    plt.plot(x, k, linewidth=2, color="green")
    for xk in zip(x, k, q):
        plt.annotate(xk[2].replace("_",""), xy=xk[:2], textcoords='data', ha='center', va="bottom")
    plt.title('' + os.path.basename(reference) + ' Ca-Ca distances', fontsize=20)
    plt.ylabel('Angstrom', fontsize=16)
    plt.xlabel('CV id', fontsize=16)
    # # plt.xlim(numpy.arange(numpy.min(x), numpy.max(x), 1.0))
    plt.xticks(numpy.arange(numpy.min(x), numpy.max(x), 1.0))
    plt.ylim((0, 4))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Alpha', markerfacecolor='r', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Beta', markerfacecolor='y', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Coil', markerfacecolor='b', markersize=10)]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.savefig(os.path.join("./", os.path.basename(reference)[:-4] + "_ca-ca_d.png"))
    plt.close()

    # r = np.arange(0, 2, 0.01)
    # theta = 2 * np.pi * r

    # ax = plt.subplot(111, projection='polar')
    v, ax = plt.subplots(1, sharex=True, figsize=(width_pic, height_pic), subplot_kw=dict(polar=True))

    # ax = axarr.subplot(111, projection='polar')
    # ax.bar(0, numpy.pi, width=0.1)
    # ax.bar(math.pi / 3.0, 3.0, width=math.pi / 3.0)

    # Adjust the axis
    # ax.set_ylim(0, numpy.pi)

    # ax.set_frame_on(False)
    # ax.axes.get_xaxis().set_visible(False)
    # ax.axes.get_yaxis().set_visible(False)

    ax.scatter([(numpy.pi / 180.0) * sole for sole in y], x, color=[cmap[u] for u in z])
    # NOTE: decomment following line to print also external angles
    ####ax.scatter([(numpy.pi/180.0 )*sole for sole in n],o, color="g")
    # NOTE: decomment following line to print also lines connecting dots
    ####ax.plot([(numpy.pi/180.0 )*sole for sole in y],x)
    # for xyq in zip([(numpy.pi/180.0 )*sole for sole in y], x, q):
    #     plt.annotate(xyq[2], xy=xyq[:2], textcoords='data')
    # ax.set_rmax(len(x))

    # ax.set_rticks([0.5, 1, 1.5, 2])  # less radial ticks
    # ax.set_rlabel_position(-22.5)  # get radial labels away from plotted line
    ax.grid(True)

    ax.set_rmax(numpy.max(x))
    legend_elements = [Line2D([0], [0], marker='o', color='w', label='Alpha', markerfacecolor='r', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Beta', markerfacecolor='y', markersize=10),
                       Line2D([0], [0], marker='o', color='w', label='Coil', markerfacecolor='b', markersize=10)]
    plt.legend(handles=legend_elements, loc='upper right')
    plt.savefig(os.path.join("./", os.path.basename(reference)[:-4] + "_Angles2.png"))
    plt.close()


#@SystemUtility.timing
def perform_superposition(reference, target, criterium_selection_core="residues", min_rmsd=0.0, max_rmsd=6.0, polyala=False, interface=False):
    """
      Superpose target pdb onto reference pdb using Characteristic Vectors annotation and graph matching.

    :param reference: The pdb reference file
    :type reference: io.TextIOWrapper
    :param target: The pdb target file
    :type target: io.TextIOWrapper
    :param min_dist_bs: Minimum distance fr CVs annotated as beta strand
    :type min_dist_bs: float
    :param max_dist_bs: Maximum distance fr CVs annotated as beta strand
    :type max_dist_bs: float
    :param strictness_ah: strictness parameter threshold for accepting ah CVs
    :type strictness_ah: float
    :param strictness_bs: strictness parameter threshold for accepting bs CVs
    :type strictness_bs: float
    :param peptide_length: Define the peptide length for computing a CV
    :type peptide_length: int
    :param write_graphml: Write graph graphml files
    :type write_graphml: bool
    :return:
    :rtype:
    """
    swap = True if len(criterium_selection_core.split("|||"))>1 and criterium_selection_core.split("|||")[1] == "swap" else False
    swapreverse = True if len(criterium_selection_core.split("|||"))>1 and criterium_selection_core.split("|||")[1] == "swapreverse" else False
    if swapreverse: swap = False
    criterium_selection_core = criterium_selection_core.split("|||")[0]

    remarks_target = None

    if isinstance(reference, tuple) and isinstance(target, tuple):
        reference, structure_reference, graph_no_coil_reference, matrix_reference, cvs_reference = reference
        target, structure_target, graph_no_coil_target, matrix_target, cvs_target = target
    else:
        structure_reference = Bioinformatics.get_structure("ref", reference)
        structure_target,remarks_target = Bioinformatics.get_structure("targ", target, get_header=True)
        # ree = reference
        # tee = target
        reference = Bioinformatics.get_list_of_residues(structure_reference, sorting=True, interface=interface)
        target = Bioinformatics.get_list_of_residues(structure_target, sorting=True, interface=interface)

        if len(reference) != len(target):
            print("SUPERPOSITION require a core of the same size: reference", len(reference), "and target", len(target), "have different size.")
            return {}

        ref_frags = [[]]
        tar_frags = [[]]

        for i in range(len(reference)):
            if i == 0 or Bioinformatics.check_continuity(reference[i], reference[i - 1]):
                ref_frags[-1].append(reference[i])
            else:
                ref_frags.append([reference[i]])

        for i in range(len(target)):
            if i == 0 or Bioinformatics.check_continuity(target[i], target[i - 1]):
                tar_frags[-1].append(target[i])
            else:
                tar_frags.append([target[i]])

        ref_frags = sorted(ref_frags, key=lambda x:len(x), reverse=True)
        tar_frags = sorted(tar_frags, key=lambda x:len(x), reverse=True)

        check_ref = [len(x) for x in ref_frags]
        check_tar = [len(x) for x in tar_frags]

        if check_ref != check_tar:
            print("SUPERPOSITION require a core of the same size: fragments in reference", check_ref,
                      "and fragments in target", check_tar, "have different sizes.")
            # tee.seek(0)
            # with open("targ" + str(time.time()) + ".pdb", "w") as gg:
            #     gg.write(tee.read())

            return {}

        map_ref = {}
        for i, ref in enumerate(ref_frags):
            if len(ref) not in map_ref:
                map_ref[len(ref)] = [i]
            else:
                map_ref[len(ref)].append(i)

        map_tar = {}
        for i,tar in enumerate(tar_frags):
            if len(tar) not in map_tar:
                map_tar[len(tar)] = [i]
            else:
                map_tar[len(tar)].append(i)

        map_ref = [ref for le,ref in map_ref.items()]
        map_tar = [tar for le,tar in map_tar.items()]

        #print(map_ref)
        #print(map_tar)
        combinations = [q for q in sorted(map_ref, key=lambda x: len(x)) for t in range(len(q))] + \
                       [q for q in sorted(map_tar, key=lambda x: len(x)) for t in range(len(q))]

        minR,mint,minu,mincore = None,None,None,None
        minrmsd = numpy.inf
        mincombi = None
        for combi in itertools.product(*combinations):
            if len(set(combi[:len(ref_frags)]))+len(set(combi[len(ref_frags):])) == len(ref_frags)*2:
                for u in range(3):
                    re1 = [(resi["N"],resi["CA"],resi["C"],resi["O"]) for cb in combi[:len(ref_frags)] for resi in ref_frags[cb][u if len(ref_frags[cb])> 5 else 0:]]
                    ta2 = [(resi["N"],resi["CA"],resi["C"],resi["O"]) for cb in combi[len(ref_frags):] for resi in tar_frags[cb][u if len(tar_frags[cb])> 5 else 0:]]
                    re1 = [atm for group in re1 for atm in group]
                    ta2 = [atm for group in ta2 for atm in group]
                    if not swap:
                        rmsd,R,t = Bioinformatics.get_rmsd_and_RT(re1, ta2, None, transform=False, n_iter=3)
                        tm = Bioinformatics.get_tm_score(re1, ta2)
                    else:
                        rmsd,R,t = Bioinformatics.get_rmsd_and_RT(ta2, re1, None, transform=False, n_iter=3)
                        tm = Bioinformatics.get_tm_score(ta2, re1)

                    ulu = u
                    core = len(re1)
                    if u > 0:
                        re1 = [(resi["N"],resi["CA"],resi["C"],resi["O"]) for cb in combi[:len(ref_frags)] for resi in ref_frags[cb][:(-u if len(ref_frags[cb])> 5 else len(ref_frags[cb]))]]
                        ta2 = [(resi["N"],resi["CA"],resi["C"],resi["O"]) for cb in combi[len(ref_frags):] for resi in tar_frags[cb][:(-u if len(tar_frags[cb])> 5 else len(ref_frags[cb]))]]
                        re1 = [atm for group in re1 for atm in group]
                        ta2 = [atm for group in ta2 for atm in group]

                        if not swap:
                            rmsd2, R2, t2 = Bioinformatics.get_rmsd_and_RT(re1, ta2, None, transform=False, n_iter=3)
                            tm2 = Bioinformatics.get_tm_score(re1, ta2)
                        else:
                            rmsd2, R2, t2 = Bioinformatics.get_rmsd_and_RT(ta2, re1, None, transform=False, n_iter=3)
                            tm2 = Bioinformatics.get_tm_score(ta2, re1)

                        #print("B",combi, rmsd2, -1*u)
                        #if combi == (0, 1, 2, 3, 4, 6, 5, 0, 1, 2, 4, 3, 6, 5):
                        #    for r in range(len(re1)):
                        #        print(r,re1[r].get_full_id(),ta2[r].get_full_id())
                        #    quit()
                        if rmsd2 < rmsd:
                            rmsd = rmsd2
                            R = R2
                            t = t2
                            ulu = -1*u
                            core = len(re1)
                            tm = tm2

                    if rmsd < minrmsd:
                        minrmsd = rmsd
                        minR = R
                        mint = t
                        minu = ulu
                        mincore = core
                        mincombi = combi
                        mintm = tm

        if minrmsd < numpy.inf:
            re1 = [(resi["N"].get_coord(),resi["CA"].get_coord(),resi["C"].get_coord(),resi["O"].get_coord()) for cb in mincombi[:len(ref_frags)] for resi in ref_frags[cb]]
            ta2 = [(resi["N"].get_coord(),resi["CA"].get_coord(),resi["C"].get_coord(),resi["O"].get_coord()) for cb in mincombi[len(ref_frags):] for resi in tar_frags[cb]]
            re1 = copy.deepcopy(numpy.array([numpy.array(atm) for group in re1 for atm in group]))
            ta2 = copy.deepcopy(numpy.array([numpy.array(atm) for group in ta2 for atm in group]))
            if not swap:
                ta2 = Bioinformatics.transform(ta2, minR, mint)
            else:
                re1 = Bioinformatics.transform(re1, minR, mint)

            # print("mincombi",mincombi)
            #print("Core rmsd was",minrmsd,end=" ")
            minrmsd = Bioinformatics.get_rmsd(re1, ta2)
            tm = Bioinformatics.get_tm_score(re1, ta2)
            #print("Final global rmsd is", minrmsd)
            # print(ref_frags)
            # print(tar_frags)
            res1 = [resi.get_full_id() for cb in mincombi[:len(ref_frags)] for resi in ref_frags[cb]]
            res2 =  [resi.get_full_id() for cb in mincombi[len(ref_frags):] for resi in tar_frags[cb]]
            # for e,r1 in enumerate(res1):
            #     print(r1,"    ",res2[e])
            # print()

        if numpy.inf > minrmsd <= minrmsd <= max_rmsd:
            print("Best RMSD achieved",minrmsd,"cutting extremities with",minu,"for a total of a core size of",mincore)
            if not swap:
                pdbtar, cnv = Bioinformatics.get_pdb_from_list_of_atoms(Bioinformatics.get_list_of_atoms(structure_target, sorting=True), renumber=False, uniqueChain=False, applyRt=(minR, mint), polyala=polyala)
            else:
                pdbtar, cnv = Bioinformatics.get_pdb_from_list_of_atoms(Bioinformatics.get_list_of_atoms(structure_reference, sorting=True), renumber=False, uniqueChain=False, applyRt=(minR, mint), polyala=polyala)
            if remarks_target is not None: pdbtar = remarks_target + "\n" + pdbtar
            return {"rmsd": minrmsd, "size": mincore, "associations": None,
                "transf": (minR,mint), "graph_ref": None, "isreversed":swapreverse,
                "grapf_targ": None, "match": None, "explored": None,
                "correlation": None, "pdb_target": pdbtar, "pdb_core_target": pdbtar,
                "ca_target": None, "tm":tm}
        elif remarks_target is not None and "REMARK ALEPH MATRIX SCORE" in remarks_target:
            scd,sca = [line for line in remarks_target.splitlines() if line.startswith("REMARK ALEPH MATRIX SCORE")][0].split()[4:6]
            scd =  float(scd)
            sca = float(sca)
            if scd >= 0.9: #<= 40:
                if not swap:
                    pdbtar, cnv = Bioinformatics.get_pdb_from_list_of_atoms(Bioinformatics.get_list_of_atoms(structure_target, sorting=True), renumber=False, uniqueChain=False, applyRt=(minR, mint), polyala=polyala)
                else:
                    pdbtar, cnv = Bioinformatics.get_pdb_from_list_of_atoms(Bioinformatics.get_list_of_atoms(structure_reference, sorting=True), renumber=False, uniqueChain=False, applyRt=(minR, mint), polyala=polyala)

                if remarks_target is not None: pdbtar = remarks_target + "\n" + pdbtar
                return {"suggested": True, "rmsd": minrmsd, "tm":tm, "size": mincore,"transf":(minR,mint), "isreversed":swapreverse, "pdb_target": pdbtar, "pdb_core_target": pdbtar, "score_mat_dist":scd, "score_mat_ang":sca}
            else:
                return {"discarded": True, "rmsd": minrmsd, "tm":tm}
        else:
            #print("RMSD Discarded",minrmsd)
            return {"discarded":True, "rmsd":minrmsd, "tm":tm}

@SystemUtility.timing
def find_local_folds_in_the_graph(reference, strictness_ah, strictness_bs, peptide_length, write_pdb, min_ah=4, min_bs=3,
                                  domains=False, compactness=False, write_graphml=False, write_graph_img=False,
                                  path_base="./"):
    # NOTE: PARAMETRIZATION=======
    # write_graphml = bool(write_graphml)
    sort_mode = "avg"

    if strictness_ah < 0:
        strictness_ah = 0.45
    if strictness_bs < 0:
        strictness_bs = 0.20

    if write_graph_img:
        path_dendo = 'os.path.basename(reference)[:-4] + "_dendo.png"'
    else:
        path_dendo = None

    weight = "distance_avg"
    pep_len = int(peptide_length)

    pdb_search_in = ""
    f = open(reference, "r")
    all_lines = f.readlines()
    f.close()
    for line in all_lines:
        if line[:6] in ["TITLE ", "CRYST1", "SCALE1", "SCALE2", "SCALE3"]:
            pdb_search_in += line

    graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
        reference, weight=weight, min_ah=min_ah, min_bs=min_bs, write_pdb=write_pdb, strictness_ah=strictness_ah,
        strictness_bs=strictness_bs, peptide_length=pep_len, domains=False, compactness=compactness)

    if domains: #To weight all the edges equally
        graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
            reference, weight=weight,  min_ah=min_ah, min_bs=min_bs, write_pdb=False, strictness_ah=strictness_ah,
            strictness_bs=strictness_bs, peptide_length=pep_len, verbose=False, domains=domains, compactness=compactness)

    graph_full_reference.vs["name"] = [str(vertex["reslist"][0][3][1]) + "_" + str(vertex["reslist"][0][2]) for vertex
                                       in graph_full_reference.vs]
    graph_no_coil_reference = graph_full_reference.vs.select(sstype_in=["ah", "bs"]).subgraph()
    eigen_ref = graph_no_coil_reference.evcent(directed=False, scale=True, weights=graph_no_coil_reference.es["weight"],
                                               return_eigenvalue=True)
    # print("EIGEN Centrality values:", eigen_ref, len(eigen_ref[0]), len(graph_no_coil_reference.vs))
    graph_no_coil_reference.vs["eigen"] = eigen_ref[0]
    graph_no_coil_reference.vs["Label"] = [d["sstype"] + "_" + d["name"] for t, d in
                                           enumerate(graph_no_coil_reference.vs)]
    #graph_no_coil_reference.write_graphml(os.path.join("./", os.path.basename(reference)[:-4] + "_graphref.graphml"))

    # NOTE: This would be the n cut chosen by the algorithm itself, is useful for community clustering decomposition but not for my application here
    ###get_community_clusters_one_step("walktrap", graph_no_coil_reference, structure_reference, pdb_search_in, os.path.basename(reference)[:-4] + "_" + str(0) + ".pdb", n=None)


    dendo = None
    vclust = None
    list_graphs = []
    dict_vclust = {}
    for cut in range(1, graph_no_coil_reference.vcount() + 1):
        if cut == 1:
            vclust, dendo = get_community_clusters_one_step("walktrap", graph_no_coil_reference, structure_reference,
                                                            pdb_search_in,
                                                            os.path.basename(reference)[:-4] + "_" + str(
                                                                cut) + ".pdb", n=cut,
                                                            print_dendo=path_dendo, return_dendo=True,
                                                            write_graph_img=write_graph_img, path_base=path_base,
                                                            write_pdb=write_pdb)
            if path_dendo:
                print("Dendogram written at:", os.path.basename(reference)[:-4] + "_dendo.png")

        else:
            vclust = get_community_clusters_one_step("walktrap", graph_no_coil_reference, structure_reference,
                                                     pdb_search_in,
                                                     os.path.basename(reference)[:-4] + "_" + str(cut) + ".pdb",
                                                     n=cut, print_dendo=None, return_dendo=False, use_dendo=dendo,
                                                     write_graph_img=write_graph_img, path_base=path_base,
                                                     write_pdb=write_pdb)
            if write_pdb:
                print("Fold extracted from dendogram cut level:", cut, "written on:",
                      os.path.basename(reference)[:-4] + "_" + str(cut) + ".pdb")

        dict_vclust[cut] = vclust
        list_graphs += vclust.subgraphs()

    if write_graphml:
        for u, g in enumerate(list_graphs):
            g_all = generate_graph_for_cvs(g, matrix_reference, cvs_reference,
                                           peptide_length=pep_len)
            g_all.write_graphml(os.path.join("./", os.path.basename(reference)[:-4] + "_" + str(u) + "_fold.graphml"))
            print("Graph written at:",
                  os.path.join("./", os.path.basename(reference)[:-4] + "_" + str(u) + "_fold.graphml"))

    return dict_vclust

def print_secondary_structure_elements(lisBigSS):
    print("=====================================================================================")
    print("      _________________________________________________________________________      ")
    for l in range(len(lisBigSS)):
        fragment = lisBigSS[l]
        print(str(fragment["pdbid"]) + "  " + str(fragment["model"]) + "  " + str(fragment["chain"]) + "  " + "[" + str(
            (fragment["resIdList"][0])[3][1]) + str((fragment["resIdList"][0])[2]) + ":" + str(
            (fragment["resIdList"][fragment["fragLength"] - 1])[3][1]) + str(
            (fragment["resIdList"][fragment["fragLength"] - 1])[2]) + "]  " + str(fragment["vecLength"]) + "  " + str(
            fragment["sstype"]) + "  " + str(fragment["sequence"]), end="")
        print(" ")
    print("      _________________________________________________________________________      ")
    print("=====================================================================================")

def __compareSeqReal(seq1, seq2, ssbridge, seq_vaso, structure):
    resids = []
    for xc in range(len(seq1)):
        if seq1[xc].upper() != "X" and seq1[xc].upper() != seq2[xc].upper():
            return False
        elif ssbridge:
            if seq1[xc].upper() == "C":
                seqv = seq_vaso[xc]
                res = Bioinformatics.get_residue(structure, seqv[1], seqv[2], seqv[3])
                resids.append(res)

    for r in range(len(resids)):
        S1 = resids[r]["SG"]
        resaS1X = float(S1.get_coord()[0])
        resaS1Y = float(S1.get_coord()[1])
        resaS1Z = float(S1.get_coord()[2])
        found = False
        for r2 in range(len(resids)):
            if r2==r: continue
            S2 = resids[r2]["SG"]

            resaS2X = float(S2.get_coord()[0])
            resaS2Y = float(S2.get_coord()[1])
            resaS2Z = float(S2.get_coord()[2])
            checkbond = numpy.sqrt(((resaS1X - resaS2X) ** 2) + ((resaS1Y - resaS2Y) ** 2) + ((resaS1Z - resaS2Z) ** 2))
            #print("-----------------------------++++++++++++++++++++++++++++=-----------------------------Checking ssbridge",checkbond,S1.get_full_id(),S2.get_full_id())
            if (2.0 <= checkbond <= 2.10):  # S-S bond is 2.05 A
                found = True
        if not found:
            return False
    # if isEqual:
    # print "Comparing",seq1,seq2
    # print "IsEqual?",isEqual
    return True

#@SystemUtility.timing
def decompose_pdb_in_cvs_mat_geom(cvs_flat, get_dot_product=False, euclidean = False, normalise = True, norm_dist=None):
    #Matrix of the mean points of a all cvs vectors
    V = numpy.array([cv[2]+((cv[3]-cv[2])/2.0) for cv in cvs_flat]) #Getting the coords of the central point in each vector a+((b-a)/2)
    if euclidean:
        V = scipy.spatial.distance.pdist(V, "euclidean")  # Getting euclidean
    else:
        V = scipy.spatial.distance.pdist(V, "sqeuclidean")   #Getting squared euclidean
    D = scipy.spatial.distance.squareform(V) #**0.5 #Getting distances

    # Normalise distances
    if norm_dist:
        D = D/norm_dist #Normalising
        D[D > 1] = 1 #Flattening values above 1

    #NOTE: Angle is given by cos(theta) = dot(a,b) / |a||b|
    # so what I can do is compute all the a and b and put in M, then compute the norm for each vectors (all the |a|) then
    # normalize each element of the vector a = (a1, a2, a3) by its norm in this way: a = (a1/|a|, a2/|a|, a3/|a|)
    # finally The Grahm matrix is the dot(M,M.T) which is already giving me the cos(theta) (just do the multiplication by hand and you will see it)

    M = numpy.array([cv[3]-cv[2] for cv in cvs_flat]) #Getting the vectors in i,j,k form
    #N = numpy.array(map(numpy.linalg.norm, M))
    N = numpy.linalg.norm(M, axis=1) #Getting the norm for each vector (a**2+b**2+c**2)**0.5
    M = M/N[:,None]  #This is magic!!! Allows me to divide every element in the ith row of matrix M by the single element of ith row in vector N
    T = numpy.dot(M,M.T) #This is computing the cos(theta)
    T[numpy.where(T < -1)] = -1.0 #fixing numerical underflow
    T[numpy.where(T > 1)]  = 1.0 #fixing numerical overflow
    Q = T.copy()
    #NOTE: Test with the ARCCOS
    if normalise:
        T = (numpy.arccos(T)*57.2958)/180.0 #This is getting theta in radians, convert it to degrees and normalize in 0-1
    else:
        T = (numpy.arccos(T)*57.2958)

    # J = N*N[:,None]
    # P = N**2
    # Q = numpy.sqrt((P*P[:,None])*(1.0-(T**2)))  # This is computing the (a**2)*(b**2)*(1-(cos_theta)**2) which the magnitude of the cross product between a and b squared
    # Q = Q/J #This should be the sin(theta)
    # print("Q")
    # print(Q)
    # print("S")
    # print(numpy.sin(numpy.arccos(T)))
    # Q = numpy.sin(numpy.arccos(T))
    # Q[numpy.where(Q < -1)] = -1.0  # fixing numerical underflow
    # Q[numpy.where(Q > 1)] = 1.0  # fixing numerical overflow
    cvs_indeces = [cv[0] for cv in cvs_flat]
    if not get_dot_product:
        return cvs_indeces,D,T #,Q
    else:
        return cvs_indeces,D,T,Q


def process_structure_fast(pdbsearch, pdbf, pattern, cvs_model, sequence, ncssearch, remove_redundance, identity, ssbridge,
                      peptide_length, connectivity, weight, strictness_ah, strictness_bs, list_chains, processed_cvlist=[], verbose=False, interface = True):

    print("Processing structure ", pdbf)
    draw_forest = False
    all_solutions = []
    strucc2 = None
    lisSS = None

    #This scaler is used to determine how many elements explore from the sorted list depending on the J threhold
    mml = MinMaxScaler(0.0, 100.0, 10000.0, 1000.0)
    max_limit = mml.scale(identity[5],integer=True)

    if not os.path.exists(pdbf):
        inpu1 = io.StringIO(SystemUtility.py2_3_unicode(pdbf))
    else:
        inpu1 = pdbf

    graph_full_reference, stru, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
        pdbsearch,
        weight=weight,
        strictness_ah=0.4,
        strictness_bs=0.3,
        peptide_length=peptide_length,
        write_pdb=False, only_reformat=True, is_model=True, verbose=False, interface=interface)


    strucc2, cvs_target = parse_pdb_and_generate_cvs(inpu1, peptide_length=peptide_length, one_model_per_nmr=True, only_reformat=not remove_redundance)
    cvs_target = [cv for cvs in cvs_target for cv in cvs]

    map_continuity = {}

    if interface: #Sorting set to False because it is unsorting when insertion code is used
        allresi = Bioinformatics.get_list_of_residues(strucc2, model=None, sorting=False, interface = True)
    else:
        allresi = Bioinformatics.get_list_of_residues(strucc2, model=None, sorting=True)
    for i in range(len(allresi)-1):
        resi1 = allresi[i].get_full_id()
        resi2 = allresi[i+1].get_full_id()
        result = Bioinformatics.check_continuity(Bioinformatics.get_residue(strucc2, resi1[1], resi1[2], resi1[3]),
                                                 Bioinformatics.get_residue(strucc2, resi2[1], resi2[2], resi2[3]))
        map_continuity[(resi1, resi2)] = result
        map_continuity[(resi2, resi1)] = result


    Fr,Dr,Tr = decompose_pdb_in_cvs_mat_geom(cvs_model)
    Ft,Dt,Tt = decompose_pdb_in_cvs_mat_geom(cvs_target)
    fra = [[g for g in group] for group in consecutive_groups(Fr)]
    fri = [[Fr.index(g) for g in group] for group in consecutive_groups(Fr)]

    #Mapping of the weights with a secondary structure meaning:
    #1) Beta strands within the same beta sheet should move less than beta strands across different sheet
    #2) Distance and angles between closest fragments

    t = (1.0/(Dr.shape[0]*Dr.shape[1]))
    #t=1.0/
    Wr = numpy.ones(Dr.shape)*t
    for frag1 in graph_full_reference.vs:
        cvs1 = sorted(frag1["cvids"])
        if len(cvs1)<=0:
            continue
        a,b = Fr.index(cvs1[0]),Fr.index(cvs1[-1]) + 1
        for frag2 in graph_full_reference.vs:
            cvs2 = sorted(frag2["cvids"])
            if len(cvs2) <= 0:
                continue

            c,d = Fr.index(cvs2[0]), Fr.index(cvs2[-1]) + 1
            if frag1["sstype"] in ["ah","bs"] and frag2["sstype"] in ["ah","bs"]:
                t = (1.0/((b-a)*(d-c)))

            if frag1["sheet"] is not None and frag2["sheet"] is not None and frag1["sheet"] == frag2["sheet"]:
                Wr[a:b,c:d] = 3*t
                Wr[a][c] = 5*t
                Wr[a][d-1] = 5*t
                Wr[b-1][c] = 5*t
                Wr[b-1][d-1] = 5*t
            # if frag1["sstype"] in ["ah","bs"] and frag2["sstype"] in ["ah","bs"] and frag1["sstype"] == frag2["sstype"]:
            #     Wr[a:b,c:d] = 0.6

    #Wr[numpy.where(numpy.abs(Tr - 1.0) <= 0.05)] += 3.0
    #Wr[numpy.where(numpy.abs(Tr + 1.0) <= 0.05)] += 3.0


    # print("Structural weights")
    #print(Wr)

    #LOOKING FOR THE SINGLE FRAGMENTS
    solutions = []
    scores = []
    #mms = MinMaxScaler(0.0, 1.0, 0.4, 1.0)

    # cutoff_d_cont = mms.scale(identity[0]/100.0)
    # cutoff_t_cont = mms.scale(identity[1]/100.0)
    # cutoff_d_jump = mms.scale(identity[4]/100.0)
    # cutoff_t_jump = mms.scale(identity[5]/100.0)
    cutoff_d_cont = identity[0] / 100.0
    cutoff_t_cont = identity[1] / 100.0
    cutoff_d_jump = identity[4] / 100.0
    cutoff_t_jump = identity[5] / 100.0
    print("Cutoff thresholds: cont dist",identity[0])
    print("Cutoff thresholds: cont angle",identity[1])
    print("Cutoff thresholds: jump dist",identity[4])
    print("Cutoff thresholds: jump angle",identity[5])
    # print("Cutoff thresholds: cont dist", identity[0], "corresponding to", cutoff_d_cont)
    # print("Cutoff thresholds: cont angle", identity[1], "corresponding to", cutoff_t_cont)
    # print("Cutoff thresholds: jump dist", identity[4], "corresponding to", cutoff_d_jump)
    # print("Cutoff thresholds: jump angle", identity[5], "corresponding to", cutoff_t_jump)

    if draw_forest: forest = igraph.Graph(directed=False)

    for i in range(len(fri)):
        bound = fri[i]
        start,end = bound[0],bound[-1]+1 #Remember to add one because the last one is not included
        ref_d = Dr[start:end,start:end]
        ref_t = Tr[start:end,start:end]
        M = ref_d.shape[0]
        ref_d = ref_d + numpy.eye(M)
        ref_t = ref_t + numpy.eye(M)

        size = end-start #Number of CVs in the fragment
        frags = [] #Fragments extracted solutions
        if verbose:
            residues = sorted(list(set([tuple(res) for cc in range(start, end) for res in cvs_reference[cc][4]])))
            print("Fragment ",i, "is", residues[0][2],residues[0][3][1],"---",residues[-1][3][1])
        for z in range(Dt.shape[0]-size+1):
            #print("start",start,"end",end,"size",size,"cut at init",z,"fin",size+z,"of",Dt.shape[0])
            s,e=z,size+z
            cmp_d = Dt[s:e,s:e]
            cmp_t = Tt[s:e,s:e]
            cmp_d = cmp_d + numpy.eye(M)
            cmp_t = cmp_t + numpy.eye(M)

            diff_d = Bioinformatics.distance_between_matrices([ref_d], [cmp_d], verbose=False)#, stat_test=False)
            diff_t = Bioinformatics.distance_between_matrices([ref_t], [cmp_t], are_cosins=True, verbose=False)#, stat_test=False)

            ####print("CONTINOUS","D",diff_d,"T",diff_t)
            if verbose:
                residues = sorted(list(set([tuple(res) for cc in range(s, e) for res in cvs_target[cc][4]])))
                print("Fragment candidate for i",i,"cvs start",s,"cvs end",e,  residues[0][2],residues[0][3],"---",residues[-1][3],
                      "score cont dist",diff_d,cutoff_d_cont,diff_d < cutoff_d_cont,
                      "score cont angle",diff_t,cutoff_t_cont,diff_t < cutoff_t_cont,
                      "result",diff_d < cutoff_d_cont or diff_t < cutoff_t_cont)

            if diff_d < cutoff_d_cont or diff_t < cutoff_t_cont: continue
            #else: print("Continue",diff_d,cutoff_d_jump,diff_t,cutoff_t_jump)
            #print("-",end="")
            frags.append([(diff_d, diff_t), (s, e)])

        frags = sorted(frags, key=lambda x: x[0], reverse=True) #Sort the frags extracted by their scores
        print("Search for fragment n. ", i, "produced",len(frags),"selected",end=" ")

        if draw_forest and i==0:
            forest.add_vertices(len(frags))
            forest.vs["nameid"] = [str(c[1][0])+"-"+str(c[1][1]) for c in frags]
            forest.vs["score_dist"] = [c[0][0] for c in frags]
            forest.vs["score_angle"] = [c[0][1] for c in frags]
            forest.vs["level"] = [i]*len(frags)

        if verbose: print("Before fragment filtering we had", len(frags))
        frags = frags[:max_limit]
        if verbose: print("After fragment filtering we have", len(frags))

        if verbose:
            print("------------------------FULL LIST OF FRAGMENTS SURVIVED----------------------")
            for frag in frags:
                print(frag)
            print("-----------------------------------------------------------------------------")

        print(len(frags),end=" ")

        if len(frags) == 0:
            break

        #TODO: Comparisons among fragments
        again = []
        agsco = []

        for lsd,frag in enumerate(frags):
            #print("Difference is", frag[0], "solu",frag[1])
            s,e=frag[1]
            current_start,current_end= fri[i][0],fri[i][-1]+1 #Remember to add one because the last one is not included
            #print("Frag. ", lsd+1, "/", len(frags))
            #print("---",current_start,current_end,"---",s,e)

            if len(solutions) > 0:
                #print("Len di solu",len(solutions))
                for solu in solutions:
                    compatible=True
                    tot_dd = 0
                    tot_tt = 0
                    ref_dd_frags = []
                    cmp_dd_frags = []
                    ref_tt_frags = []
                    cmp_tt_frags = []
                    ww_frags = []
                    for n,pos in enumerate(solu):
                        g,f=pos
                        #print("Evaluating",pos)
                        start,end=fri[n][0],fri[n][-1]+1
                        #print("g,f",g,f,"start,end",start,end,"current_start,current_end",current_start,current_end)
                        if verbose:
                            residues_a = sorted(list(set([tuple(res) for cc in range(s, e) for res in cvs_target[cc][4]])))
                            residues_b = sorted(list(set([tuple(res) for cc in range(g, f) for res in cvs_target[cc][4]])))
                            print("----------------------Compatibility of a fragment with a previous solution---------------------------")
                            print("Frag ",lsd, residues_a[0][2],residues_a[0][3],"---",residues_a[-1][3],"with Solu",n,residues_b[0][2],residues_b[0][3],"---",residues_b[-1][3])
                            # print("Frag ",lsd, residues_a[0][2],residues_a[0][3][1],"---",residues_a[-1][3][1],"with Solu",n,residues_b[0][2],residues_b[0][3][1],"---",residues_b[-1][3][1])
                        if len(set(list(range(g, f))) & set(list(range(s, e)))) > 0:
                            if verbose:
                                print("Not acceptable because cvs in common:",set(list(range(g,f)))&set(list(range(s,e))))
                                print("----------------------------------------------------------------------------------")

                            compatible = False
                            break

                        if connectivity:
                            info_g = (cvs_target[g][4][0][2], cvs_target[g][4][0][3][1])
                            info_s = (cvs_target[s][4][0][2], cvs_target[s][4][0][3][1])
                            info_start = (cvs_target[start][4][0][2], cvs_target[start][4][0][3][1])
                            info_current_start = (cvs_target[current_start][4][0][2], cvs_target[current_start][4][0][3][1])
                            #print("g,s",info_g, info_s,"start,current_start",info_start,info_current_start)

                            if info_g[0] == info_s[0] and info_start[0] == info_current_start[0]:
                                if numpy.sign(info_g[1] - info_s[1]) != numpy.sign(info_start[1] - info_current_start[1]):
                                    if verbose:
                                        print("Incompatible because connectivity test 1, it is not respected: ",numpy.sign(info_g[1] - info_s[1]), numpy.sign(info_start[1] - info_current_start[1]))
                                        print(
                                            "----------------------------------------------------------------------------------")

                                    compatible = False
                                    break
                            elif numpy.sign(ord(info_g[0])-ord(info_s[0])) != numpy.sign(ord(info_start[0])-ord(info_current_start[0])):
                                if verbose:
                                    print("Incompatible because connectivity test 2, it is not respected: ",
                                                  numpy.sign(ord(info_g[0])-ord(info_s[0])),
                                                  numpy.sign(ord(info_start[0])-ord(info_current_start[0])))
                                    print(
                                        "----------------------------------------------------------------------------------")

                                compatible = False
                                break

                        if interface:
                            sss = sorted(list(set([tuple(res) for cc in list(range(g, f)) + list(range(s, e)) for res in cvs_target[cc][4]])),  key=lambda x: (x[2], x[3][2], x[3][1]))
                            ttt = sorted(list(set([tuple(res) for cc in list(range(start, end)) + list(range(current_start, current_end)) for res in cvs_model[cc][4]])), key=lambda x: (x[2], x[3][2], x[3][1]))
                        else:
                            sss = sorted(list(set([tuple(res) for cc in sorted(list(range(g, f)) + list(range(s, e))) for res in cvs_target[cc][4]])), key=lambda x: (x[2], x[3][1]))
                            ttt = sorted(list(set([tuple(res) for cc in sorted(list(range(start, end)) + list(range(current_start, current_end))) for res in cvs_model[cc][4]])), key=lambda x: (x[2], x[3][1]))

                        if len(sss) != len(ttt):
                            if verbose:
                                print("Incompatible because have different number of residues",len(sss),len(ttt))
                                print("----------------------------------------------------------------------------------")

                            compatible = False
                            break

                        tar_frags = [[]]
                        for v in range(len(sss)):
                            # print((sss[v][:-1], sss[v - 1][:-1]))
                            if v == 0 or ((sss[v][:-1], sss[v - 1][:-1]) in map_continuity and map_continuity[(sss[v][:-1], sss[v - 1][:-1])]):
                                # print('esta')
                                tar_frags[-1].append(sss[v])
                            else:
                                # print('no esta')
                                tar_frags.append([sss[v]])

                        ref_frags = sorted([fri[n][-1] - fri[n][0] + 3, fri[i][-1] - fri[i][0] + 3], reverse=True)
                        tar_frags = sorted([len(p) for p in tar_frags], reverse=True)

                        # print("CCC")
                        # print("ref",ref_frags)
                        # print("tar",tar_frags)
                        # print("compatible",compatible)
                        if ref_frags != tar_frags:
                            if verbose:
                                print("Incompatible because the reference sizes are not compatible with the target sizes",ref_frags,tar_frags)
                                print("----------------------------------------------------------------------------------")

                            compatible = False
                            break
                        # print("DDD",compatible)

                        cmp_dd = Dt[g:f,s:e]
                        cmp_tt = Tt[g:f,s:e]
                        ref_dd = Dr[start:end,current_start:current_end]
                        ref_tt = Tr[start:end,current_start:current_end]
                        ref_ww = Wr[start:end,current_start:current_end]
                        ref_dd_frags.append(ref_dd.copy())
                        #print("ref",ref_dd_frags)
                        ref_tt_frags.append(ref_tt.copy())
                        cmp_dd_frags.append(cmp_dd.copy())
                        #print("cmp",cmp_dd_frags)
                        cmp_tt_frags.append(cmp_tt.copy())
                        ww_frags.append(ref_ww.copy())
                        #print("N",n,"Ref",(g,f,s,e))
                        #print("N",n,"Cmp",(start,end,current_start,current_end))

                        # fre1_dd = Dr[start:end, start:end]
                        # fre1_tt = Tr[start:end, start:end]
                        # fre2_dd = Dr[current_start:current_end, current_start:current_end]
                        # fre2_tt = Tr[current_start:current_end, current_start:current_end]
                        # fcm1_dd = Dt[g:f, g:f]
                        # fcm1_tt = Tt[g:f, g:f]
                        # fcm2_dd = Dt[s:e, s:e]
                        # fcm2_tt = Tt[s:e, s:e]
                        # fre1_dd = fre1_dd + numpy.eye(fre1_dd.shape[0])
                        # fre1_tt = fre1_tt + numpy.eye(fre1_tt.shape[0])
                        # fre2_dd = fre2_dd + numpy.eye(fre2_dd.shape[0])
                        # fre2_tt = fre2_tt + numpy.eye(fre2_tt.shape[0])
                        # fcm1_dd = fcm1_dd + numpy.eye(fcm1_dd.shape[0])
                        # fcm1_tt = fcm1_tt + numpy.eye(fcm1_tt.shape[0])
                        # fcm2_dd = fcm2_dd + numpy.eye(fcm2_dd.shape[0])
                        # fcm2_tt = fcm2_tt + numpy.eye(fcm2_tt.shape[0])

                    if compatible:
                        #print("I",i,len(ref_dd_frags))
                        fixed_t = Bioinformatics.distance_between_matrices(ref_tt_frags, cmp_tt_frags, weights_list=ww_frags, are_cosins=True, verbose=False)#, stat_test=True if i==2 and n==1 else False, verbose=True if i==2 and n==1 else False)
                        fixed_d = Bioinformatics.distance_between_matrices(ref_dd_frags, cmp_dd_frags, weights_list=ww_frags, verbose=False)#, stat_test=True if i==2 and n==1 else False, verbose=True if i==2 and n==1 else False)

                        # fixed_t = Bioinformatics3.distance_between_matrices([ref_tt],[cmp_tt], weights_list=[ref_ww], are_cosins=True, verbose= False)#, stat_test=True if i==2 and n==1 else False, verbose=True if i==2 and n==1 else False)
                        # fixed_d = Bioinformatics3.distance_between_matrices([ref_dd],[cmp_dd], weights_list=[ref_ww], verbose= False)#, stat_test=True if i==2 and n==1 else False, verbose=True if i==2 and n==1 else False)

                        # conti1_d = Bioinformatics3.distance_between_matrices(fre1_dd,fcm1_dd)
                        # conti1_t = Bioinformatics3.distance_between_matrices(fre1_tt,fcm1_tt, are_cosins=True)
                        # conti2_d = Bioinformatics3.distance_between_matrices(fre2_dd,fcm2_dd)
                        # conti2_t = Bioinformatics3.distance_between_matrices(fre2_tt,fcm2_tt, are_cosins=True)
                        tot_dd = fixed_d
                        tot_tt = fixed_t
                        if fixed_d < cutoff_d_jump or fixed_t < cutoff_t_jump:
                            if verbose:
                                print("Incompatible because fixed fragment threhsold are not repsected distance:",fixed_d,cutoff_d_jump,"angles",fixed_t,cutoff_t_jump)
                                print("----------------------------------------------------------------------------------")
                            compatible=False
                            #NOTE comment this break if we move the indent out of this for
                            ###break
                        #else:
                            #print("True fixed",fixed_d,cutoff_d_jump,"fixedt",fixed_t,cutoff_t_jump)

                        # print("BBB")
                        # if i == 2 and n==1:
                        # print("\n","N",n,"JUMPS D", fixed_d, "T", fixed_t)
                        # print(set(list(range(g,f)))&set(list(range(s,e))))
                        # print(sss)
                        # print(ttt)
                        # print("DISTANCES")
                        # Bioinformatics3.distance_between_matrices(ref_dd, cmp_dd, are_cosins=False, verbose=True, weights=Wr)
                        # print("ANGLES")
                        # Bioinformatics3.distance_between_matrices(ref_tt, cmp_tt, are_cosins=True, verbose=True, weights=Wr)

                        ####print("\nJUMPS D", fixed_d, "T", fixed_t, "CONTINOUS D1", conti1_d, "T1", conti1_t, "D2", conti2_d, "T2", conti2_t, "TOT D", tot_dd, "TOT T", tot_tt)

                    ####NOTE: HEre it was the ident
                    #print("DDD")
                    #print("RRR",compatible)
                    if compatible:
                        # for ll in alld:
                        #     print(ll)
                        # print("----",tot_dd/len(solu))
                        #print("compatible",solu)
                        o = copy.deepcopy(solu)
                        o.append((s,e))
                        again.append(((tot_dd, tot_tt),o))
                        if verbose:
                            print("The combination has been accepted!!!!")
                            print("----------------------------------------------------------------------------------")

                        if draw_forest:
                            named = "_".join([str(jj[0])+"-"+str(jj[1]) for jj in solu])
                            for node in forest.vs.select(level=0):
                                z = [p for p in forest.get_all_shortest_paths(node.index) if len(p)==len(solu)]
                                names = ["_".join([forest.vs[t]["nameid"] for t in c]) for c in z]
                                #print(names)
                                if named in names:
                                    q = names.index(named)
                                    #print(names,"found at",q,names[q],z[q])
                                    last = z[q][-1]
                                    vid = forest.vcount()
                                    forest.add_vertex(nameid=str(s)+"-"+str(e),score_dist=tot_dd,score_angle=tot_tt,level=i)
                                    forest.add_edge(vid,last)

                        #print(".",end="")
            else:
                again.append((frag[0],[frag[1]]))
        #print("")

        if len(again) == 0:
            print("but not a valid combination considering fixed fragments")
            break

        again = sorted(again, key=lambda x: x[0], reverse=True)[:max_limit*2]
        if i==len(fri)-1: #and len(again) > max_limit: #i>1:
            exten = [sorted([ran for fr in solu[1] for ran in range(fr[0],fr[1])]) for solu in again]
            correct = [exten[0]]
            lich = [0]
            for ex in range(1,len(exten)):
                #print("Evaluating",ex,"/",len(exten))
                insert = True
                for px in correct:
                    sss = set([tuple(res) for cc in exten[ex] for res in cvs_target[cc][4]])
                    ppp = set([tuple(res) for cc in px for res in cvs_target[cc][4]])

                    if len(sss-ppp) < (5*len(sss))/100.0:
                        if verbose:
                            print("--------")
                            print("The current solution is sharing more than the 5 per cent with another solution so it is removed as considered redundant")
                            print("SET 1",sss,"SET 2",ppp)
                            print("--------")
                        insert = False
                        break
                if insert:
                    correct.append(exten[ex])
                    lich.append(ex)
        else:
            lich = range(len(again))

        if verbose: print("Before solution filtering we had ", len([d for p,d in enumerate(again) if p in lich]))
        again = [d for p,d in enumerate(again) if p in lich][:max_limit]
        if verbose: print("After solution filtering we have", len(again))

        if verbose:
            print("-------------------FULL LIST OF SOLUTIONS SURVIVED--------------------")
            for ag in again:
                print(ag)
            print("----------------------------------------------------------------------")

        #
        # if i == 5:
        #        for ag in again:
        #            print(ag)
        dictio=[]
        agsco = [ag[0] for ag in again] #This line is fundamental that is done before the next one
        again = [ag[1] for ag in again]

        #print("++++++++++++++++++++++++++++++++++++")
        #print(again)


        solutions = again
        scores = agsco
        print("and number of valid combinations considered fixed fragments:",len(solutions))

    all_solutions = []

    chains_dict={}

    for tr, solu in enumerate(solutions):
        chains_dict[tr] = {}
        if len(solu) == len(fri):
            for num, wert in enumerate(solu):
                for cc in range(wert[0], wert[1]):
                    for res in cvs_target[cc][4]:
                        if len(list_chains) == len(fri):
                            chains_dict[tr][Bioinformatics.get_residue(strucc2, res[1], res[2], res[3]).get_unpacked_list()[0].get_full_id()[:-1]] = list_chains[num]
                        else:
                            chains_dict[tr][Bioinformatics.get_residue(strucc2, res[1], res[2], res[3]).get_unpacked_list()[0].get_full_id()[:-1]] = 'A'

    for tr,solu in enumerate(solutions):
        if len(solu) == len(fri):
            residues = sorted(list(set([tuple(res) for wert in solu for cc in range(wert[0],wert[1]) for res in cvs_target[cc][4]])))
            cvs_p = [cvs_target[cc] for wert in solu for cc in range(wert[0],wert[1])]
            back_atoms_list = []
            atoms_list = []
            chains_list = []

            if len(sequence) > 0:
                seqr = "".join([Bioinformatics.AADICMAP[Bioinformatics.get_residue(strucc2, resi[1], resi[2], resi[3]).get_resname().upper()] for resi in residues])
                # print(sequence)
                # print(seqr)
                # print()
                if not __compareSeqReal(sequence, seqr, ssbridge, residues, strucc2):
                    if verbose: print("Solution ",residues[0][2], residues[0][3][1], "---",residues[-1][3][1], "with sequence", seqr, "do not align according sequence schema", sequence)
                    continue

            for resi in residues:  back_atoms_list += Bioinformatics.get_backbone(Bioinformatics.get_residue(strucc2, resi[1], resi[2], resi[3]), without_CB=False)
            for resi in residues:  atoms_list += Bioinformatics.get_residue(strucc2, resi[1], resi[2], resi[3]).get_unpacked_list()

            ida = atoms_list[0].get_full_id()
            all_solutions.append((back_atoms_list, atoms_list, cvs_p, ida[0], ida[1], scores[tr]))

    #all_solutions = all_solutions[:1]
    if draw_forest:
        forest.write_graphml(os.path.join("./", os.path.basename(pdbf)[:-4] + "_forest.graphml"))

    return strucc2, all_solutions, chains_dict, []


def start_new_process_pdb(pdbsearch, pdbf, idname, DicParameters, doCluster=True, getStructure=False, justBackbone=True, processed_cvlist=[], size_tree=5, interface = False, verbose=False):
    try:
        nombre = "undefined"
        toremove = ""

        if not os.path.exists(pdbsearch):
            inpu1 = io.StringIO(SystemUtility.py2_3_unicode(pdbsearch))
        else:
            inpu1 = pdbsearch

        if not os.path.exists(pdbf):
            inpu2 = io.StringIO(SystemUtility.py2_3_unicode(pdbf))
            nombre = idname
        else:
            inpu2 = pdbf
            nombre = os.path.basename(inpu2)[:-4]

        if idname == "undefined":
            idname = nombre

        if DicParameters["exclude_sequence"]:
            structure_read = Bioinformatics.get_structure("pdb", pdbf)
            seqbychain = {chai.get_id(): "".join(
                [Bioinformatics.AADICMAP[d.get_parent().get_resname()] for d in Bioinformatics.get_list_of_atoms(structure_read, type_atom="CA") if
                 d.get_parent().get_resname() in Bioinformatics.AADICMAP and d.get_parent().get_full_id()[2] == chai.get_id()]) for
                 chai in Bio.PDB.Selection.unfold_entities(structure_read, "C")}

            chain_to_skip = set([])
            for cid, seqc in seqbychain.items():
                z = galign(seqc,DicParameters["exclude_sequence"])
                q = sum([1.0 for g in range(len(z[0][0])) if z[0][0][g] != "-" and z[0][1][g] != "-" and z[0][0][g]==z[0][1][g]])/len(z[0][0]) if len(z) > 0 else 0.0
                if q >= 0.9:
                    print("Skipping chain",cid,"Because global alignment with exclude_sequence is",q)
                    chain_to_skip.add(cid)

            if len(chain_to_skip) == len(seqbychain.keys()):
                print("All the chains have been skipped so this pdb is abandoned!")
                if getStructure:
                    return ("", [], [], None, [])
                else:
                    return ("", [], [], [])
            elif len(chain_to_skip) > 0:
                print("Pdb have been reduced to only chains: ",set(list(seqbychain.keys()))-chain_to_skip)
                atoms = [atm for atm in Bioinformatics.get_list_of_atoms(structure_read) if atm.get_full_id()[2] not in chain_to_skip]
                pdbftext,_ = Bioinformatics.get_pdb_from_list_of_atoms(atoms)
                with tempfile.NamedTemporaryFile(mode="w", delete=False) as tff:
                    tff.write(pdbftext)
                    toremove = tff.name
                    pdbf = tff.name
            else:
                print("No chains have been removed from the analysis.")

        print("================= The size of the fragment for computing a CV is:", DicParameters["peptide_length"], "aa.======================")
        list_chains = []
        for fragment in DicParameters['listFragments']:
            if fragment['sstype']== 'ah' or fragment['sstype'] =='bs':
                list_chains.append(fragment['reslist'][0][2])
        list_chains = sorted(list_chains)

        structure, solutions, chains_dict, lisBigSS = process_structure_fast(pdbsearch, pdbf, DicParameters["pattern"], DicParameters["cvsModel"],
                                                           DicParameters["sequenceLocate"], DicParameters["ncssearch"],
                                                           DicParameters["remove_redundance"],
                                                           DicParameters["identity"],
                                                           DicParameters["ssbridge"], DicParameters["peptide_length"],
                                                           DicParameters["connectivity"], DicParameters["weight"],
                                                           DicParameters["strictness_ah"], DicParameters["strictness_bs"],list_chains,
                                                           processed_cvlist=processed_cvlist, interface=interface, verbose=verbose)
        for p in range(len(solutions)):
            # NOTE: TEMPORANEO
            # Is this "for" block  really needed
            pda = solutions[p][0]



        # if len(solutions) > 0:
        #     dic_json = {pdbf.split('/')[-1]: len(solutions)}
        #     modify_json_file('output.json', dic_json, lib_generation_extraction=True)

        start = int(math.floor(numpy.sqrt(len(solutions) / 2.0)))
        #TODO: When finally rewrite clusterization and introduce Kmean through sklearn end should be end=2*start
        end = start #2*start
        #TODO: Substitute cluster_by_kmeans by a simple function that is only preparing the necessary things to continue further without this kmeans clustering
        clusts = prepare_for_grouping(DicParameters, "", pdbsearch, solutions, DicParameters["strictness_ah"], DicParameters["strictness_bs"], DicParameters["peptide_length"],
                                      structure=structure, writeOutput=False, criterium_selection_core=DicParameters["criterium_selection_core"])

        newSolList = []
        newSolRed = []
        for tuplos in clusts:
            solu_chain = ''
            if len(tuplos) == 0:
                continue
            cosa = tuplos[0]
            name = cosa[:3]
            pdbid, model, IdSolution = name
            #print("choosing: ", pdbid, model, IdSolution)
            q = 0
            remark = "REMARK ALEPH MATRIX SCORE "+str(solutions[IdSolution][5][0])+" "+str(solutions[IdSolution][5][1])+"\n"

            if justBackbone:
                ppp = Bioinformatics.get_pdb_from_list_of_atoms(solutions[IdSolution][0], renumber=True, dictio_chains = chains_dict[IdSolution], interface=interface)
                # ppp = Bioinformatics.get_pdb_from_list_of_atoms(solutions[IdSolution][0], renumber=False,dictio_chains=chains_dict[IdSolution]) #For not renaming the residues
                # ppp = Bioinformatics.get_pdb_from_list_of_atoms(solutions[IdSolution][0], renumber=False) #For not renaming the residues and the chains
                for atom_i, atom in enumerate(solutions[IdSolution][0]):
                    if atom_i == 0:
                        res_chain = atom.get_full_id()[2]
                        res_icode = atom.get_full_id()[3][2]
                    elif atom.get_full_id()[2] != res_chain or atom.get_full_id()[3][2] != res_icode:
                        solu_chain = 'interface'
                        break
                    else:
                        solu_chain = 'monomer'

            else:
                ppp = Bioinformatics.get_pdb_from_list_of_atoms(solutions[IdSolution][1], renumber=True, dictio_chains = chains_dict[IdSolution], polyala=False)
                # ppp = Bioinformatics3.get_pdb_from_list_of_atoms(solutions[IdSolution][1], renumber=False,dictio_chains=chains_dict[IdSolution]) #For not renaming the residues
                # ppp = Bioinformatics.get_pdb_from_list_of_atoms(solutions[IdSolution][0], renumber=False) #For not renaming the residues and the chains
                for atom_i, atom in enumerate(solutions[IdSolution][1]):
                    if atom_i == 0:
                        res_chain = atom.get_full_id()[2]
                        res_icode = atom.get_full_id()[3][2]
                    elif atom.get_full_id()[2] != res_chain or atom.get_full_id()[3][2] != res_icode:
                        solu_chain = 'interface'
                        break
                    else:
                        solu_chain = 'monomer'
                q = 1
            ppp = (remark+ppp[0],ppp[1], solu_chain)
            newSolList.append(ppp)
            newSolRed.append(map(lambda x: Bioinformatics.get_pdb_from_list_of_atoms(solutions[x[:3][2]][q]), tuplos[1:]))

        if os.path.exists(toremove): os.remove(toremove)

        if getStructure:
            return (idname, newSolList, lisBigSS, structure, newSolRed)
        else:
            return (idname, newSolList, lisBigSS, newSolRed)
    except:
        print("An error occured while parsing or decoding PDB ", idname, " PID: " + str(os.getpid()))
        print(sys.exc_info())
        traceback.print_exc(file=sys.stdout)
        t = datetime.datetime.now()
        epcSec = time.mktime(t.timetuple())
        now = datetime.datetime.fromtimestamp(epcSec)
        print("" + str(now.ctime()) + "\tError parsing or decoding: " + str(idname) + " PID: " + str(os.getpid()) + "\n" + str(sys.exc_info()) + "\n")
        # quit()

        if os.path.exists(toremove): os.remove(toremove)

        if getStructure:
            return ("", [], [], None, [])
        else:
            return ("", [], [], [])

def create_pdbs(wdir, pdbsearch, pdbfstru, solutions, solred, pdbid, model, DicParameters, representative=False, superpose=True, superpose_exclude=1, return_pdbstring=False, nilges=10, ref_pdbid = False,
                superposition_mode=False, pdb_model=None, interface=False,  justBackbone=True):

    global number_of_solutions

    if not return_pdbstring and not os.path.exists(os.path.join(wdir, "./library/")):
        os.makedirs(os.path.join(wdir, "./library"))
    # if not return_pdbstring:
    #     if os.path.exists(os.path.join(wdir, "./library/")):
    #         shutil.rmtree("./library/")
    #     os.makedirs(os.path.join(wdir, "./library"))

    allpdbfstru = None
    allpdbtar = None
    if isinstance(pdbfstru, str) and os.path.exists(pdbfstru):
        f = open(pdbfstru, "r")
        allpdbtar = f.read()
        allpdbfstru = allpdbtar.splitlines()
        f.close()
    else:
        allpdbtar = pdbfstru.read()
        allpdbfstru = allpdbtar.splitlines()

    title = ""
    for linea in allpdbfstru:
        lis = linea.split()
        if len(lis) > 1 and lis[0] == "TITLE":
            title += " ".join(lis[1:]) + "  "

    pdbsol = []
    allred = {}
    total = len(solutions)
    dictio_pdb_best = {}

    if superposition_mode:
        pdbid = os.path.basename(pdb_model)[:-17]

    for p in range(len(solutions)):
        # NOTE: TEMPORANEO
        # Is this "for" block  really needed
        pda = solutions[p][0]
        # pda = ""
        # for lineaA in solutions[p][0].splitlines():
        #     if lineaA.startswith("ATOM") or lineaA.startswith("HETATM"):
        #         for lineaB in allpdbfstru:
        #             if lineaA[30:54] in lineaB and lineaA[16:20] in lineaB:
        #                 pda += lineaB
        #     elif lineaA.startswith("REMARK"):
        #         pda += lineaA
        # print pda
        write = True
        # NOTE: Activate next for debugging
        # superpose = False
        listrmsd = []
        if superpose:
            write = False
            # NOTE: Substituted the old getSuperimp with the new ALEPH.perform_superposition
            # {"rmsd": best_rmsd, "size": best_size, "associations": asso, "transf": best_transf, "graph_ref": g_a, "grapf_targ": g_b, "correlation": corr,"pdb_target":pdbmod}
            # TODO: inside the functions that opens reference and target it must be checked if it is a path or not to not open it
            dictio_super = perform_superposition(reference=io.StringIO(SystemUtility.py2_3_unicode(pdbsearch)),
                                                 target=io.StringIO(SystemUtility.py2_3_unicode(pda)),
                                                 criterium_selection_core=DicParameters["criterium_selection_core"],
                                                 min_rmsd=DicParameters["minRMSD"], max_rmsd=DicParameters["maxRMSD"], interface=interface)

            name = "" + str(pdbid) + "_" + str(model) + "_" + str(p) + ".pdb"
            if "discarded" in dictio_super and dictio_super["discarded"]:
                ###print("Superposition discarded because rmsd",dictio_super["rmsd"],"is greater than threshold.")
                rmsT = dictio_super["rmsd"]
                name = "" + str(pdbid) + "_" + str(model) + "_" + str(p) + ".pdb"
                pda = solutions[p][0]
                write = False
            elif "suggested" in dictio_super and dictio_super["suggested"]:
                rmsT = dictio_super["rmsd"]
                if not dictio_super["isreversed"]:
                    name = "" + str(pdbid) + "_" + str(model) + "_" + str(p) + ".pdb"
                else:
                    name = "" + "output" + "_" + str(model) + "_" + str(p) + ".pdb"

                pda = dictio_super["pdb_target"]
                size = dictio_super["size"]
                scd = dictio_super["score_mat_dist"]
                sca = dictio_super["score_mat_ang"]
                if not os.path.exists(os.path.join(wdir, "./library/suggested")): os.makedirs(os.path.join(wdir, "./library/suggested"))
                if not dictio_super["isreversed"]:
                    with open(os.path.join(wdir, "./library/suggested/" + name), "w") as j:
                        j.write(pda)
                else:
                    pdacc, cnv = Bioinformatics.get_pdb_from_list_of_atoms(Bioinformatics.get_list_of_atoms(Bioinformatics.get_structure("result", io.StringIO(SystemUtility.py2_3_unicode(allpdbtar))), sorting=True), renumber=False, uniqueChain=False, applyRt=dictio_super["transf"])
                    rmsT = dictio_super["rmsd"]

                print("SUPERPOSITION", name, rmsT, dictio_super["tm"],"stored in ./library/suggested because score matrix distance is", scd)
                write = False
            elif not dictio_super:
                ###print("Cannot superpose", name)
                rmsT = 100
                name = "" + str(pdbid) + "_" + str(model) + "_" + str(p) + ".pdb"
                pda = solutions[p][0]
                write = False
            else:
                if not dictio_super["isreversed"]:
                    pdacc = dictio_super["pdb_target"]
                    rmsT = dictio_super["rmsd"]
                else:
                    pdacc, cnv = Bioinformatics.get_pdb_from_list_of_atoms(Bioinformatics.get_list_of_atoms(Bioinformatics.get_structure("result", io.StringIO(SystemUtility.py2_3_unicode(allpdbtar))), sorting=True), renumber=False, uniqueChain=False, applyRt=dictio_super["transf"], polyala= justBackbone)
                    rmsT = dictio_super["rmsd"]
                    name = "" + "output" + "_" + str(model) + "_" + str(p) + ".pdb"

            if DicParameters["minRMSD"] <= rmsT <= DicParameters["maxRMSD"]:
                if representative and pdbid in dictio_pdb_best and rmsT < dictio_pdb_best[pdbid][1]:
                    os.remove(dictio_pdb_best[pdbid][0])
                    dictio_pdb_best[pdbid] = (os.path.join(wdir, "./library/" + name), rmsT)
                    print("SUPERPOSITION", name, rmsT, dictio_super["tm"], dictio_super["size"], title)
                    listrmsd.append((rmsT, name, dictio_super["size"], dictio_super["tm"]))
                    write = True
                elif representative and pdbid not in dictio_pdb_best:
                    dictio_pdb_best[pdbid] = (os.path.join(wdir, "./library/" + name), rmsT)
                    print("SUPERPOSITION", name, rmsT, dictio_super["tm"], dictio_super["size"], title)
                    listrmsd.append((rmsT, name, dictio_super["size"], dictio_super["tm"]))
                    write = True
                elif not representative:
                    print("SUPERPOSITION",name, rmsT, dictio_super["tm"], dictio_super["size"], title)
                    listrmsd.append((rmsT, name, dictio_super["size"], dictio_super["tm"]))
                    write = True

        if write:
            dic_json = {}
            if not return_pdbstring:
                if not dictio_super["isreversed"]:
                    name = "" + str(pdbid) + "_" + str(model) + "_" + str(p) + ".pdb"
                else:
                    name = "" + "output" + "_" + str(model) + "_" + str(p) + ".pdb"

                f = open(os.path.join(wdir, "./library/" + name), "w")
                f.write('REMARK   ' + solutions[p][2] + '\n' + pdacc) #pda
                f.close()
                number_of_solutions += 1
                f = open(os.path.join(wdir, "library/" + "list_rmsd.txt"), "a")
                listrmsd = sorted(listrmsd)
                for pair in listrmsd:
                    f.write(pair[1] + "  " + str(pair[0]) + "  " + str(pair[3]) + "  " + str(pair[2]) + "\n")
                    #NOTE: For Ana: in this way you are overwriting dic_json at each iteration is that what you want?
                    dic_json = {pair[1]: {'rmsd and core': str(pair[0]) + "  " + str(pair[2]),'tmscore':str(pair[3])}}
                f.close()
                #modify_json_file('output.json', dic_json, lib_generation_superposition=True)
            else:
                if not os.path.isdir(os.path.join(wdir, "library")):
                    os.makedirs(os.path.join(wdir, "library"))
                f = open(os.path.join(wdir, "library/" + "list_rmsd.txt"), "a")
                listrmsd = sorted(listrmsd)
                for pair in listrmsd:
                    f.write(pair[1] + "  " + str(pair[0]) + "  " + str(pair[3]) + "  " + str(pair[2]) + "\n")
                    dic_json = {pair[1]: {'rmsd and core': str(pair[0]) + "  " + str(pair[2]),'tmscore':str(pair[3])}}
                f.close()
                #modify_json_file('output.json', dic_json, lib_generation_superposition=True)

                if not dictio_super["isreversed"]:
                    name = "" + str(pdbid) + "_" + str(model) + "_" + str(p) + ".pdb"
                else:
                    name = "" + "output" + "_" + str(model) + "_" + str(p) + ".pdb"
                pdbsol.append((name, pdacc)) #pda
                allred["" + str(pdbid) + "_" + str(model) + "_" + str(p)] = []
                for pda2 in solred[p]:
                    if not dictio_super["isreversed"]:
                        name2 = "" + str(pdbid) + "_" + str(model) + "_" + str(total) + ".pdb"
                    else:
                        name2 = "" + "output" + "_" + str(model) + "_" + str(total) + ".pdb"
                    allred["" + str(pdbid) + "_" + str(model) + "_" + str(p)].append((name2, pda2))
                    total += 1

    return pdbsol, allred

def start_new_process(wdir, pdbsearch, pdbstruc, DicParameters, doCluster, backbone, cvs_list_str, superpose, idn,
                      return_pdbstring, thresh, superpose_exclude, superposition_mode, pdb_model, interface = False, verbose=False):
    # NOTE: Activate next line for debugging
    # doCluster = False
    try:
        niceness = os.nice(0)
        os.nice(5 - niceness)
    except:
        pass


    if idn is None or len(idn) == 0: idn="undefined"

    (idname, solList, lisBigSS, solred) = start_new_process_pdb(pdbsearch, pdbstruc, idn, DicParameters, doCluster=doCluster,
                                                                justBackbone=backbone, processed_cvlist=cvs_list_str, interface=interface, verbose=verbose)

    if len(lisBigSS) == 0:
        lisBigSS = [{"model": 0}]

    if len(solList) > 0:
        if len(cvs_list_str) > 0:
            idname = idn
            if return_pdbstring:
                return create_pdbs(wdir, pdbsearch, io.StringIO(SystemUtility.py2_3_unicode(pdbstruc)), solList, solred, idname,
                                   lisBigSS[0]["model"], DicParameters,
                                   superpose=superpose, return_pdbstring=return_pdbstring,
                                   superpose_exclude=superpose_exclude, nilges=DicParameters["nilges"],
                                   representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)

            else:
                if doCluster:
                    pdbsols, solred = create_pdbs(wdir, pdbsearch, io.StringIO(SystemUtility.py2_3_unicode(pdbstruc)), solList, solred, idname,
                                                  lisBigSS[0]["model"], DicParameters,
                                                  superpose=superpose, return_pdbstring=True,
                                                  superpose_exclude=superpose_exclude, nilges=DicParameters["nilges"],
                                                  representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)
                    cluster_by_rmsd(os.path.join(wdir, "./library"), pdbsols, solred, thresh, superpose_exclude, DicParameters["nilges"], DicParameters["weight"],
                                    DicParameters["strictness_ah"], DicParameters["strictness_bs"], DicParameters["peptide_length"],
                                    criterium_selection_core=DicParameters["criterium_selection_core"])
                else:
                    create_pdbs(wdir, pdbsearch, io.StringIO(SystemUtility.py2_3_unicode(pdbstruc)), solList, solred, idname,
                                lisBigSS[0]["model"], DicParameters,
                                superpose=superpose, return_pdbstring=False, superpose_exclude=superpose_exclude,
                                nilges=DicParameters["nilges"], representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbon)

        else:
            if os.path.exists(pdbstruc):
                if return_pdbstring:
                    return create_pdbs(wdir, pdbsearch, pdbstruc, solList, solred, idname, lisBigSS[0]["model"],
                                       DicParameters, superpose=superpose,
                                       return_pdbstring=return_pdbstring, superpose_exclude=superpose_exclude,
                                       nilges=DicParameters["nilges"], representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)

                else:
                    if doCluster:
                        pdbsols, solred = create_pdbs(wdir, pdbsearch, pdbstruc, solList, solred, idname,
                                                      lisBigSS[0]["model"], DicParameters,
                                                      superpose=superpose, return_pdbstring=True,
                                                      superpose_exclude=superpose_exclude, nilges=DicParameters["nilges"],
                                                      representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)

                        cluster_by_rmsd(os.path.join(wdir, "./library"), pdbsols, solred, thresh, superpose_exclude,
                                        DicParameters["nilges"], DicParameters["weight"],
                                        DicParameters["strictness_ah"], DicParameters["strictness_bs"],
                                        DicParameters["peptide_length"],
                                        criterium_selection_core=DicParameters["criterium_selection_core"])
                    else:
                        create_pdbs(wdir, pdbsearch, pdbstruc, solList, solred, idname, lisBigSS[0]["model"],
                                    DicParameters, superpose=superpose,
                                    return_pdbstring=False, superpose_exclude=superpose_exclude, nilges=DicParameters["nilges"],
                                    representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)


            else:
                idname = idn
                if return_pdbstring:
                    return create_pdbs(wdir, pdbsearch, io.StringIO(SystemUtility.py2_3_unicode(pdbstruc)), solList, solred, idname,
                                       lisBigSS[0]["model"], DicParameters,
                                       superpose=superpose, return_pdbstring=return_pdbstring,
                                       superpose_exclude=superpose_exclude, nilges=DicParameters["nilges"],
                                       representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)

                else:
                    if doCluster:
                        pdbsols, solred = create_pdbs(wdir, pdbsearch, io.StringIO(SystemUtility.py2_3_unicode(pdbstruc)), solList, solred,
                                                      idname, lisBigSS[0]["model"],
                                                      DicParameters, superpose=superpose, return_pdbstring=True,
                                                      superpose_exclude=superpose_exclude, nilges=DicParameters["nilges"],
                                                      representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)

                        cluster_by_rmsd(os.path.join(wdir, "./library"), pdbsols, solred, thresh, superpose_exclude,
                                        DicParameters["nilges"], DicParameters["weight"],
                                        DicParameters["strictness_ah"], DicParameters["strictness_bs"],
                                        DicParameters["peptide_length"],
                                        criterium_selection_core=DicParameters["criterium_selection_core"])
                    else:
                        create_pdbs(wdir, pdbsearch, io.StringIO(SystemUtility.py2_3_unicode(pdbstruc)), solList, solred, idname,
                                    lisBigSS[0]["model"], DicParameters,
                                    superpose=superpose, return_pdbstring=False, superpose_exclude=superpose_exclude,
                                    nilges=DicParameters["nilges"], representative=DicParameters["representative"], superposition_mode=superposition_mode, pdb_model=pdb_model, interface=interface, justBackbone=backbone)

    else:
        return [], {}

def evaluate_pdb(sym, wdir, pdbf, cvs_list_str, pdbstruc, strucc, lisBig, i, pattern, pattern_cvs, highd,
                 doCluster, superpose, process_join, pdbsearch, pdbn, thresh, superpose_exclude,
                 peptide_length, sequence, ncssearch, multimer, c_angle, c_dist,
                 c_angle_dist, c_cvl_diff, j_angle, j_dist, j_angle_dist, j_cvl_diff, rmsd_min, rmsd_max,
                 ssbridge, connectivity, nilges, enhance_fold, representative, pdb_model,
                 sidechains, weight, strictness_ah, strictness_bs,
                 criterium_selection_core,  exclude_sequence, superposition_mode, return_pdbstring=False, interface=False, verbose=False):
    global number_of_solutions

    excluded_sequence = None
    if exclude_sequence:
        excluded_sequence = ""
        with open(exclude_sequence, "r") as f:
            for line in f.readlines():
                if not line.startswith(">"): excluded_sequence += line[:-1]

    #lisBig = sorted(lisBig, key=lambda x: x["fragLength"], reverse=True)
    DicParameters = {}
    DicParameters["nameExecution"] = "WOW"
    DicParameters["structure"] = strucc
    DicParameters["cvsModel"] = pattern_cvs
    DicParameters["pattern"] = pattern
    DicParameters["listFragments"] = lisBig
    DicParameters["highest_distance"] = highd
    DicParameters["sequenceLocate"] = sequence
    DicParameters["exclude_sequence"] = excluded_sequence
    DicParameters["ncssearch"] = ncssearch
    DicParameters["remove_redundance"] = multimer
    DicParameters["identity"] = (
        float(c_cvl_diff), float(c_angle), float(c_dist), float(c_angle_dist),  float(j_cvl_diff), float(j_angle),
        float(j_dist), float(j_angle_dist))
    DicParameters["minRMSD"] = float(rmsd_min)
    DicParameters["maxRMSD"] = float(rmsd_max)
    DicParameters["ssbridge"] = bool(ssbridge)
    DicParameters["connectivity"] = bool(connectivity)
    DicParameters["nilges"] = int(nilges)
    DicParameters["enhance_fold"] = bool(enhance_fold)
    DicParameters["peptide_length"] = int(peptide_length)
    DicParameters["representative"] = bool(representative)
    DicParameters["weight"] = weight
    DicParameters["strictness_ah"] = strictness_ah
    DicParameters["strictness_bs"] = strictness_bs
    DicParameters["criterium_selection_core"] = criterium_selection_core

    # Case .pdb file

    if os.path.exists(pdbf):
        print("Evaluating structure n.", i + 1, pdbf)
    else:
        print("Evaluating structure n.", i + 1)

    if not return_pdbstring:
        # Case real .pdb process the matrix
        if cvs_list_str is None and pdbstruc is None:
            if pdbn == "":
                p = sym.spawn_function_with_multiprocessing(target=start_new_process, args=(
                    wdir, pdbsearch, pdbf, DicParameters, doCluster, not sidechains, [],
                    superpose, "", False, thresh, superpose_exclude, superposition_mode, pdb_model,  interface, verbose))
            else:
                p = sym.spawn_function_with_multiprocessing(target=start_new_process, args=(
                    wdir, pdbsearch, pdbf, DicParameters, doCluster, not sidechains, [],
                    superpose, pdbn, False, thresh, superpose_exclude, superposition_mode, pdb_model, interface, verbose))

        if process_join:
            p.join()
    else:
        if cvs_list_str is None and pdbstruc is None:
            if pdbn == "":
                return start_new_process(wdir, pdbsearch, pdbf, DicParameters, doCluster,
                                         not sidechains, [], superpose, "", True, thresh,
                                         superpose_exclude, superposition_mode, pdb_model, interface, verbose)
            else:
                return start_new_process(wdir, pdbsearch, pdbf, DicParameters, doCluster,
                                         not sidechains, [], superpose, pdbn, True, thresh,
                                         superpose_exclude, superposition_mode, pdb_model, interface, verbose)

def evaluate_model(pdbmodel, enhance_fold, peptide_length, weight, strictness_ah, strictness_bs):
    graph_full_reference, strucc, pattern, pattern_cvs, highd = annotate_pdb_model_with_aleph(pdbmodel,
                                                                                              weight=weight,
                                                                                              strictness_ah=strictness_ah,
                                                                                              strictness_bs=strictness_bs,
                                                                                              peptide_length=peptide_length,
                                                                                              write_pdb=False, is_model=True)

    all_frags = get_all_fragments(graph_full_reference)

    pdbsearch = Bioinformatics.get_pdb_from_list_of_frags("0", all_frags, strucc, "", externalRes=[], normalize=False)[1]
    full_list = []
    lisBig = sorted(all_frags, key=lambda x: x["fragLength"], reverse=True)

    print_secondary_structure_elements(lisBig)

    if enhance_fold:
        if lisBig[-1]["fragLength"] % 2 == 0:
            peptide_length = lisBig[-1]["fragLength"] - 1
        else:
            peptide_length = lisBig[-1]["fragLength"] - 2

        graph_full_reference, strucc, pattern, pattern_cvs, highd = annotate_pdb_model_with_aleph(pdbmodel,
                                                                                                  weight=weight,
                                                                                                  strictness_ah=strictness_ah,
                                                                                                  strictness_bs=strictness_bs,
                                                                                                  peptide_length=peptide_length,
                                                                                                  write_pdb=False, is_model=True)
        all_frags = get_all_fragments(graph_full_reference)

        pdbsearch = Bioinformatics.get_pdb_from_list_of_frags("0", all_frags, strucc, "", externalRes=[], normalize=False)[1]
        full_list = []
        lisBig = sorted(all_frags, key=lambda x: x["fragLength"], reverse=True)

        print_secondary_structure_elements(lisBig)

    print_pattern(pattern)

    return pattern, pattern_cvs, highd, pdbsearch, strucc, lisBig, peptide_length

@SystemUtility.timing
@SystemUtility.deprecated('Must be changed soon for a better cluster algorithm not based in references')
def cluster_by_rmsd(directory, pdbsol, solred, thresh, superpose_exclude, nilges, weight, strictness_ah, strictness_bs, peptide_length, criterium_selection_core="residues", min_rmsd=0.0, max_rmsd=6.0):
    howmany = 0
    if not os.path.exists(directory):
        os.makedirs(directory)

    print("Start to clustering...", len(pdbsol))
    print("redundancy: ", len(solred))
    for root, subFolders, files in os.walk(directory):
        for fileu in files:
            if fileu.startswith("rmsd_list"):
                howmany += 1
                list_ref = []
                f = open(os.path.join(root, fileu), "r")
                allfiles = f.readlines()
                f.close()
                reference = None
                for line in allfiles:
                    if line.startswith("Reference"):
                        reference = os.path.join(root, line.split()[1])
                        list_ref.append(reference)

                #NOTE: before list of fragments was used to guide superpos,
                # graph_full_reference, stru, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(reference,
                #                                                                                             weight=weight,
                #                                                                                             strictness_ah=strictness_ah,
                #                                                                                             strictness_bs=strictness_bs,
                #                                                                                             peptide_length=peptide_length,
                #                                                                                             write_pdb=False)
                # all_frags = get_all_fragments(graph_full_reference)
                for reference in list_ref:
                    dds = open(reference, "r")
                    modelpdbstring = dds.read()
                    dds.close()
                    remainings = []
                    for it in range(len(pdbsol)):
                        item = pdbsol[it]
                        name, pda = item
                        listona = name.split("_")
                        pdbid = listona[0][:4]
                        model = listona[1]
                        IdSolution = listona[2]
                        if "." in list(IdSolution):
                            IdSolution, ext = IdSolution.split(".")
                        nomefile = os.path.join(root, str(pdbid) + "_" + str(model) + "_" + str(IdSolution) + ".pdb")
                        #NOTE: Substituted the old getSuperimp with the new ALEPH.perform_superposition
                        # {"rmsd": best_rmsd, "size": best_size, "associations": asso, "transf": best_transf, "graph_ref": g_a, "grapf_targ": g_b, "correlation": corr,"pdb_target":pdbmod}
                        #TODO: inside the functions that opens reference and target it must be checked if it is a path or not to not open it
                        dictio_super = perform_superposition(reference=io.StringIO(SystemUtility.py2_3_unicode(modelpdbstring)), target=io.StringIO(SystemUtility.py2_3_unicode(pda)),
                                                       criterium_selection_core=criterium_selection_core, min_rmsd=0.0, max_rmsd=100)

                        if not dictio_super:
                            print("Cannot superpose", name)
                            rmsdVALFar = 100
                        else:
                            pdacc = dictio_super["pdb_target"]
                            rmsdVALFar = dictio_super["rmsd"]
                            rmsdVALFar = dictio_super["rmsd"]
                            print(rmsdVALFar, str(pdbid) + "_" + str(model) + "_" + str(IdSolution) + ".pdb", '\t', reference)


                        # print "---",pdbid,model,IdSolution,rmsdVALFar
                        if rmsdVALFar <= thresh:
                            flo = open(os.path.join(root, fileu), "a")
                            flo.write(str(pdbid) + "_" + str(model) + "_" + str(IdSolution) + "\t" + str(rmsdVALFar) + "\t" + "Reference: " + os.path.basename(reference) + "\n")
                            flo.close()
                            fla = open(nomefile, "w")
                            fla.write(pdacc)
                            fla.close()
                            # print "Search redundant",""+str(pdbid)+" "+str(model)+" "+str(IdSolution)
                            for item2 in solred["" + str(pdbid) + "_" + str(model) + "_" + str(IdSolution)]:
                                name2, pda2 = item2
                                listona2 = name2.split("_")
                                pdbid2 = listona2[0][:4]
                                model2 = listona2[1]
                                IdSolution2 = listona2[2]
                                if "." in list(IdSolution2):
                                    IdSolution2, ext2 = IdSolution2.split(".")
                                nomefile2 = os.path.join(root, str(pdbid2) + "_" + str(model2) + "_" + str(IdSolution2) + ".pdb")
                                flo2 = open(os.path.join(root, fileu), "a")
                                flo2.write(str(pdbid2) + " " + str(model2) + " " + str(IdSolution2) + " --\n")
                                flo2.close()
                                fla2 = open(nomefile2, "w")
                                fla2.write(pda2[0])
                                fla2.close()
                        else:
                            remainings.append((name, pdacc))

                    pdbsol = remainings
    howmany += 1
    if len(pdbsol) > 0:
        while len(pdbsol) > 0:
            remainings = []
            listmodel = []
            modelpdbstring = ""
            for i in range(len(pdbsol)):
                item = pdbsol[i]
                name, pda = item
                listona = name.split("_")
                pdbid = listona[0][:4]
                model = listona[1]
                IdSolution = listona[2]
                #print(listona)
                if "." in list(IdSolution):
                    IdSolution, ext = IdSolution.split(".")

                if i == 0:
                    nameref = name
                    writePath = os.path.join(directory, str(howmany))
                    if not os.path.exists(writePath):
                        os.makedirs(writePath)

                    howmany += 1
                    flo = open(os.path.join(writePath, "rmsd_list"), "a")
                    flo.write("Reference: " + os.path.basename(name) + "\n")
                    flo.close()
                    fla = open(os.path.join(writePath, name), "w")
                    fla.write(pda)
                    fla.close()

                    # NOTE: before list of fragments was used to guide superpos,
                    # graph_full_reference, stru, matrix_reference, cvs_reference, highd_reference= annotate_pdb_model_with_aleph(os.path.join(writePath, name),
                    #                                                                                             weight=weight,
                    #                                                                                             strictness_ah=strictness_ah,
                    #                                                                                             strictness_bs=strictness_bs,
                    #                                                                                             peptide_length=peptide_length,
                    #                                                                                             write_pdb=False)
                    # all_frags = get_all_fragments(graph_full_reference)
                    modelpdbstring = pda
                    # print "Reference",os.path.basename(name)
                    # print "Search redundant ",""+str(pdbid)+" "+str(model)+" "+str(IdSolution)
                    for item2 in solred["" + str(pdbid) + "_" + str(model) + "_" + str(IdSolution)]:
                        name2, pda2 = item2
                        listona2 = name2.split("_")
                        pdbid2 = listona2[0][:4]
                        model2 = listona2[1]
                        IdSolution2 = listona2[2]
                        if "." in list(IdSolution2):
                            IdSolution2, ext2 = IdSolution2.split(".")
                        nomefile2 = os.path.join(writePath,
                                                 str(pdbid2) + "_" + str(model2) + "_" + str(IdSolution2) + ".pdb")
                        # if rmsdVALFar2 <= thresh:
                        # print "write: ",str(pdbid2)+" "+str(model2)+" "+str(IdSolution2)
                        flo2 = open(os.path.join(writePath, "rmsd_list"), "a")
                        flo2.write(str(pdbid2) + "_" + str(model2) + "_" + str(IdSolution2) + " --\n")
                        flo2.close()
                        fla2 = open(nomefile2, "w")
                        fla2.write(pda2[0])
                        fla2.close()
                    continue

                nomefile = os.path.join(writePath, str(pdbid) + "_" + str(model) + "_" + str(IdSolution) + ".pdb")
                # {"rmsd": best_rmsd, "size": best_size, "associations": asso, "transf": best_transf, "graph_ref": g_a, "grapf_targ": g_b, "correlation": corr,"pdb_target":pdbmod}
                dictio_super = perform_superposition(reference=io.StringIO(SystemUtility.py2_3_unicode(modelpdbstring)),
                                                     target=io.StringIO(SystemUtility.py2_3_unicode(pda)),
                                                     criterium_selection_core=criterium_selection_core,
                                                     min_rmsd=0.0, max_rmsd=100)
                if not dictio_super:
                    print("Cannot superpose", name)
                    rmsdVALFar = 100
                else:
                    pdacc = dictio_super["pdb_target"]
                    rmsdVALFar = dictio_super["rmsd"]

                if rmsdVALFar <= thresh:
                    # print pdbid,model,IdSolution,rmsdVALFar,thresh
                    flo = open(os.path.join(writePath, "rmsd_list"), "a")
                    flo.write(str(pdbid) + "_" + str(model) + "_" + str(IdSolution) + "\t" + str(rmsdVALFar) + "\t" + "Reference: " + nameref + "\n")
                    #flo.write(str(pdbid) + "_" + str(model) + "_" + str(IdSolution) + "\t" + str(rmsdVALFar) + "\t" + str(os.path.basename(name)) + "\n")
                    flo.close()
                    fla = open(nomefile, "w")
                    fla.write(pdacc)
                    fla.close()
                    # print "Search redundant",""+str(pdbid)+" "+str(model)+" "+str(IdSolution)
                    for item2 in solred["" + str(pdbid) + "_" + str(model) + "_" + str(IdSolution)]:
                        name2, pda2 = item2
                        listona2 = name2.split("_")
                        pdbid2 = listona2[0][:4]
                        model2 = listona2[1]
                        IdSolution2 = listona2[2]
                        if "." in list(IdSolution2):
                            IdSolution2, ext2 = IdSolution2.split(".")
                        nomefile2 = os.path.join(writePath,
                                                 str(pdbid2) + "_" + str(model2) + "_" + str(IdSolution2) + ".pdb")

                        # if rmsdVALFar2 <= thresh:
                        # print "write: ",str(pdbid2)+" "+str(model2)+" "+str(IdSolution2)
                        flo2 = open(os.path.join(writePath, "rmsd_list"), "a")
                        flo2.write(str(pdbid2) + " " + str(model2) + " " + str(IdSolution2) + " --\n")
                        flo2.close()
                        fla2 = open(nomefile2, "w")
                        fla2.write(pda2[0])
                        fla2.close()

                else:
                    remainings.append((name, pda))

            pdbsol = remainings

def cluster_by_rmsd_range(wdir, directory_database, n_clusters ,n_ranges):
    list_rmsd = []
    print("\n============Starting clustering by rmsd_range algorithm============")
    if not os.path.isfile(os.path.join(wdir, "library/list_rmsd.txt")):
        print("There are no models contained in the library. Ending.")
        return
    else:
        f = open(os.path.join(wdir, directory_database, "list_rmsd.txt"), "r")
        lines = f.readlines()
        f.close()
        for line in lines:
            line = (line.strip()).split()
            list_rmsd.append((line[0], float(line[1])))
        list_rmsd = sorted(list_rmsd, key=lambda pdb: pdb[1])
        #print(list_rmsd)
        if len(lines) < 10: return

        if n_ranges >= len(lines):
            print("The number of files extracted is less than the number of requested ranges. Number of ranges set to:", len(lines))
            n_ranges = len(lines)
        if n_clusters >= len(lines):
            print("The number of files extracted is less than the number of requested clusters. Number of clusters set to:", len(lines))
            n_clusters = len(lines)
        clusters = chunkIt(list_rmsd, n_ranges)

        select_random_pdb_from_list(clusters, directory_database, n_clusters, n_ranges, wdir)

def chunkIt(seq, num):
    # Check valores float???
    if len(seq) < num:
        print("Ther are more number of ranges than files, number of ranges set to: ", len(seq))
        num = len(seq)
    avg = len(seq) / float(num)
    out = []
    last = 0.0

    while last < len(seq):
        out.append(seq[int(last):int(last + avg)])
        last += avg

    return out

def select_random_pdb_from_list(list_pdb, directory_database, n_cluster, n_ranges, wdir):
    included = 0  # Number of models included in the library
    cluster = n_ranges - 1
    total_sup = 0 #Number of pdbs which have passed the superposition filtering

    if n_ranges > 1:
        for element in list_pdb:
            total_sup += len(element)
    else:
        total_sup = len(list_pdb)


    if not os.path.exists('library_cluster'):
        os.makedirs('library_cluster')

    print("n_cluster",n_cluster, "n_ranges",n_ranges)
    while included < n_cluster:
        if n_ranges > 1:
            select_pdb = random.sample(list_pdb[cluster], 1)
            select_pdb = select_pdb[0]
            list_pdb[cluster].remove(select_pdb)
            if not os.path.exists(os.path.join(wdir, directory_database, select_pdb[0])): continue
            print("Selected:", os.path.join(wdir, directory_database, select_pdb[0]))

            if cluster > 0:
                cluster -= 1
            else:
                cluster = n_ranges - 1
        else:
            select_pdb = random.sample(list_pdb, 1)
            list_pdb.remove(select_pdb[0])
            if not os.path.exists(os.path.join(wdir, directory_database, select_pdb[0])): continue
            print("Selected:", os.path.join(wdir, directory_database, select_pdb[0]))

        shutil.move(os.path.join(wdir, directory_database, select_pdb[0]),
                    os.path.join(wdir, "library_cluster", select_pdb[0]))
        included += 1

    dic_json = {'pdb superposed': total_sup,'number clusters': n_cluster}
    modify_json_file('output.json', dic_json, clustering=True)

    shutil.copy(os.path.join(wdir, directory_database, 'list_rmsd.txt'),
                os.path.join(wdir, "library_cluster", 'list_rmsd.txt'))

    tar_dir(os.path.join(wdir, directory_database), os.path.join(wdir,  os.path.basename(directory_database)+'_lib.tar.gz'))

def tar_dir(dir_path, tar_name):
    with tarfile.open(tar_name, "w:gz") as tar_handle:
        for root, dirs, files in os.walk(dir_path):
            for file in files:
                tar_handle.add(os.path.join(root, file))
    shutil.rmtree(dir_path, ignore_errors=True)

@SystemUtility.timing
def prepare_for_grouping(DicParameters, directory, referenceFile, heapSolutions, strictness_ah, strictness_bs, peptide_length,
                         structure=None,
                         writeOutput=True, criterium_selection_core="residues"):
    import decimal
    decimal.getcontext().prec = 2

    nameExp = DicParameters["nameExecution"]

    data = []
    ref_labels = []
    if heapSolutions is None:
        n_errors = 0
        nsoluzioni = 1
        for root, subFolders, files in os.walk(directory):
            for fileu in files:
                pdbf = os.path.join(root, fileu)
                if pdbf.endswith(".pdb"):
                    print("Evaluating file ", nsoluzioni, " that is ", pdbf)
                    nodo = os.path.basename(pdbf)[:-4]
                    pdbid = nodo.split("_")[0]
                    model = nodo.split("_")[1]
                    IdSolution = nodo.split("_")[2]
                    error = False
                    rmsd = 0.0
                    # {"rmsd": best_rmsd, "size": best_size, "associations": asso, "transf": best_transf, "graph_ref": g_a, "grapf_targ": g_b, "correlation": corr,"pdb_target":pdbmod}
                    dictio_super = perform_superposition(reference=referenceFile, target=pdbf,
                                                         criterium_selection_core=criterium_selection_core,
                                                         min_rmsd=0.0, max_rmsd=100)
                    if not dictio_super:
                        n_errors += 1
                        rmsd = -100
                        error = True
                    else:
                        rmsd = dictio_super["rmsd"]
                        if rmsd < 0:
                            n_errors += 1
                            error = True

                    if error:
                        print("Found an rmsd error for: ", pdbf)
                    nsoluzioni += 1
                    structure = Bioinformatics.get_structure("" + str(pdbid) + "_" + str(model) + "_" + str(IdSolution), pdbf)
                    data.append([rmsd])
                    ref_labels.append((pdbid, model, IdSolution))
                    print("processed node", pdbid, model, IdSolution, rmsd, "n_sol:", nsoluzioni, "n_err:", n_errors)
    elif heapSolutions is not None and structure is not None:
        cvs_model = DicParameters["cvsModel"]

        listaSolutions = []
        IdSolution = 0
        back_atm_li = None
        rmsT = None
        pdbid = None
        model = None

        for solu in heapSolutions:
            back_atm_li, atm_li, cvs, pdbid, model, score = solu
            data.append([])
            ref_labels.append((pdbid, model, IdSolution))
            # print "Preparing",pdbid,model,IdSolution,rmsT,tensorInertia,com,shape_par
            IdSolution += 1

    subclu = []
    if len(data) > 0:
        subclu = [[] for _ in range(len(data))]
        for p in range(len(data)):
            subclu[p].append((ref_labels[p]) + tuple(data[p]))

        if writeOutput:
            f = open("./clusters.txt", "w")
            f.write("Number of clusters: " + str(len(subclu)) + "\n\n")
            for clus in subclu:
                f.write("===================== " + str(len(clus)) + "\n")
                for pdbin in clus:
                    f.write("\t")
                    f.write("NAME: " + str(pdbin[0]) + " " + str(pdbin[1]) + " " + str(pdbin[2]) + " ")
                    f.write("RMSD: " + str(pdbin[3]) + " ")
                    f.write("Tensor of inertia: " + ("%.2f" % pdbin[4]) + " " + ("%.2f" % pdbin[5]) + " " + (
                        "%.2f" % pdbin[6]) + " ")
                    f.write("Center of mass: " + ("%.2f" % pdbin[7]) + " " + ("%.2f" % pdbin[8]) + " " + (
                        "%.2f" % pdbin[9]) + " ")
                    f.write("Shape Par.: " + ("%.2f" % pdbin[10]) + " " + ("%.2f" % pdbin[11]) + " " + (
                        "%.2f" % pdbin[12]) + " ")
                    f.write("\n")
                f.write("=====================\n")
            f.close()

    nodisalvati = 1
    # print "len(subclu)",len(subclu)
    return subclu

#@SystemUtility.timing
def evaluate_grid_job(idname):
    global toclientdir
    global number_of_solutions
    global doCluster_global

    print("Evaluating job ", idname)
    os.remove(os.path.join(toclientdir, idname + ".tar.gz"))
    p = subprocess.Popen('grep "_0_" ' + os.path.join(toclientdir, idname + ".out"), shell=True, stdin=subprocess.PIPE,
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    outp, errp = p.communicate()
    outp = outp.decode("ascii")
    errp = errp.decode("ascii")
    outp = outp.strip()
    f = open(os.path.join(toclientdir, "rmsd_list.txt"), "a")
    f.write(outp + "\n")
    f.close()
    os.remove(os.path.join(toclientdir, idname + ".out"))

    while 1:
        try:
            # members = tar.getmembers()
            # names = tar.getnames()
            # print "members",members
            # print "names",names
            tar = tarfile.open(os.path.join(toclientdir, idname + "_res.tar.gz"), "r:gz")
            infile = tar.extractfile(idname + "_res_out.data")
            break
        except:
            continue

    # nPDBs = pickle.load(infile)
    thresh = pickle.load(infile)
    superpose_exclude = pickle.load(infile)
    nilges = pickle.load(infile)
    weight = pickle.load(infile)
    strictness_ah = pickle.load(infile)
    strictness_bs = pickle.load(infile)
    peptide_length =  pickle.load(infile)
    criterium_selection_core = pickle.load(infile)

    library = os.path.join(toclientdir, "../library/")
    if not os.path.exists(library):
        os.makedirs(library)

    collpdbs = pickle.load(infile)
    solred = pickle.load(infile)
    infile.close()
    tar.close()

    if doCluster_global:
        cluster_by_rmsd(library, collpdbs, solred, thresh, superpose_exclude, nilges, weight, strictness_ah, strictness_bs, peptide_length, criterium_selection_core=criterium_selection_core)
    else:
        for item in collpdbs:
            pathp = os.path.join(library, item[0])
            f = open(pathp, "w")
            f.write(item[1])
            f.close()

    try:
        os.remove(os.path.join(toclientdir, idname + "_res_out.data"))
    except:
        pass
    try:
        os.remove(os.path.join(toclientdir, idname + "_res.tar.gz"))
    except:
        pass
    try:
        os.remove(os.path.join(toclientdir, idname + ".sh"))
    except:
        pass

    return True


def prepare_and_launch_job(cm, baseline, idname, pdbsfiles, supercomputer, pdb_model):
    global toclientdir
    global PATH_NEW_BORGES
    global PATH_PYTHON_INTERPRETER
    global NUMBER_OF_PARALLEL_GRID_JOBS

    # print "Number of jobs to evaluate:  "+str(len(listjobs))

    pdbf = os.path.join(toclientdir, idname + ".tar")
    fro = open(os.path.join(toclientdir, idname + "PARAM"), "wb")
    pickle.dump(len(pdbsfiles), fro)
    for valo in pdbsfiles:
        pdbf_n, pdbf_c, cvs_list_str, pdbstruc = valo
        pickle.dump(pdbf_n, fro)
        pickle.dump(pdbf_c, fro)
        pickle.dump(cvs_list_str, fro)
        pickle.dump(pdbstruc, fro)
    fro.close()
    tar = tarfile.open(pdbf, "a")
    tar.add(os.path.join(toclientdir, idname + "PARAM"), arcname=idname + "PARAM")
    tar.close()
    os.remove(os.path.join(toclientdir, idname + "PARAM"))
    compri = gzip.open(pdbf + ".gz", 'wb')
    fion = open(pdbf, "rb")
    compri.write(fion.read())
    fion.close()
    compri.close()
    os.remove(pdbf)

    if supercomputer is not None and os.path.exists(supercomputer):
        comando = "nice -n5 " + PATH_PYTHON_INTERPRETER+ " " + PATH_NEW_BORGES + " -j " + os.path.join(toclientdir,
                                                                        str(idname) + ".tar.gz") + " " + " ".join(
            baseline) + " > " + os.path.join(toclientdir, str(idname) + ".out")
        SystemUtility.launchCommand(comando, os.path.join(toclientdir, str(idname) + ".out"),
                                    """Job ended with success""", 1, single=True)
    else:
        print("# of jobs to evaluate: ", len(SystemUtility.LISTJOBS))
        while len(SystemUtility.LISTJOBS) > NUMBER_OF_PARALLEL_GRID_JOBS:
            time.sleep(3)

        if hasattr(cm, "channel"):
            cm.copy_local_file(pdbf + ".gz", pdbf + ".gz", force_cumulative=False)

        job = Grid.gridJob(str(idname))

        if hasattr(cm, "channel"):
            job.setInitialDir(cm.get_remote_pwd())
        else:
            job.setInitialDir(os.path.abspath(toclientdir))

        script = """#!/bin/bash
if [ -f """ + PATH_NEW_BORGES + """ ]; then
        """ + PATH_PYTHON_INTERPRETER+ " " + PATH_NEW_BORGES + """ """ + """ """.join(baseline) + """ --targz """ + str(idname) + """.tar.gz
else
        """ + PATH_PYTHON_INTERPRETER+ " " + os.path.basename(PATH_NEW_BORGES) + """ """ + """ """.join(baseline) + """ --targz """ + str(idname) + """.tar.gz
fi
"""
        f = open(os.path.join(toclientdir, str(idname) + ".sh"), "w")
        f.write(script)
        f.close()

        st = os.stat(os.path.join(toclientdir, str(idname) + ".sh"))
        os.chmod(os.path.join(toclientdir, str(idname) + ".sh"),
                 st.st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

        job.setExecutable(os.path.join(toclientdir, str(idname) + ".sh"))
        job.addInputFile(str(idname) + ".tar.gz", False)
        job.addInputFile(os.path.basename(pdb_model), False)
        #job.addInputFile(PATH_NEW_BORGES, False)
        job.addOutputFile(str(idname) + ".out", False)
        # job.addOutputFile(str(idname)+"_res.tar.gz",False)
        (nc, nq) = cm.submitJob(job, isthelast=True)
        SystemUtility.LISTJOBS[idname] = [(os.path.join(toclientdir, str(idname) + ".out"), "success", 1, True, "")]
        print("Search in " + str(nq) + " structure submitted to the cluster " + str(nc))


@SystemUtility.timing
def cluster_library(directory_database="./library", rmsd_clustering=1.5, clustering_mode='rmsd_range', number_of_ranges=7000,
                    number_of_clusters=7000, exclude_residues_superpose=3, ssbridge=None, nilges=None, work_directory="./"):

    wdir = work_directory
    if clustering_mode == 'rmsd_range':
        number_of_clusters = int(number_of_clusters)
        number_of_ranges = int(number_of_ranges)
        while len(multiprocessing.active_children()) > 0: pass
        cluster_by_rmsd_range(wdir, directory_database, number_of_clusters, number_of_ranges)
    elif clustering_mode == 'random_sampling':
        number_of_clusters = int(number_of_clusters)
        while len(multiprocessing.active_children()) > 0: pass
        print("\n============Starting clustering by random_sampling algorithm============")
        if not os.path.isfile(os.path.join(wdir, os.path.join(directory_database,"list_rmsd.txt"))):
            print("There are no models contained in the library. Ending.")
        else:
            pdb_files = []
            for root, subFolders, files in os.walk(os.path.join(wdir, directory_database)):
                for fileu in files:
                    pdbf = os.path.join(root, fileu)
                    if pdbf.endswith(".pdb"):
                        pdb_files.append(os.path.basename(pdbf))
            if number_of_clusters > len(pdb_files):
                number_of_clusters = len(pdb_files)
                print(
                    "The number of files extracted is less than the number of requested clusters. Number of clusters set to:",
                    len(pdb_files))
            select_random_pdb_from_list(pdb_files, directory_database, number_of_clusters, 1, wdir)
    elif clustering_mode == 'rmsd':
        #TODO: I should put here the clustering code by rmsd and doing that in a similar way to ex compact_library so have a look at that
        while len(multiprocessing.active_children()) > 0: pass
        print("\n============Starting writing clustered_library by rmsd============")
        list_ref = []
        total_sup = 0  # Number of pdbs that have passed the filtering check in the superposition
        total_ref = 0  # Number of pdb that are reference of a cluster (== number of clusters)
        if not os.path.isdir(os.path.join(wdir, directory_database)):
            print("There are no models contained in the library. Ending.")
        else:
            g = open(wdir + 'list_rmsd_clustering.txt', 'a')
            for root, subFolders, files in os.walk(os.path.join(wdir, directory_database)):
                for fileu in files:
                    if fileu.startswith("rmsd_list"):
                        # howmany += 1
                        # list_ref = []
                        f = open(os.path.join(root, fileu), "r")
                        allfiles = f.readlines()
                        f.close()
                        # reference = None
                        for line in allfiles:
                            total_sup += 1
                            g.write(line)
                            if line.startswith("Reference"):
                                total_ref += 1
                                reference = os.path.join(root, line.split()[1])
                                list_ref.append(reference)
                        os.remove(os.path.join(root, fileu))
            g.close()
            try:
                dic_json = {'pdb superposed': total_sup,
                            'number clusters': total_ref}
                modify_json_file('output.json', dic_json, clustering=True)
            except:
                print("Fail in writing json file")

            if not os.path.exists('library_cluster'):
                os.makedirs('library_cluster')
            shutil.move(os.path.join(wdir, 'list_rmsd_clustering.txt'),
                        os.path.join(wdir, "library_cluster", 'list_rmsd_clustering.txt'))

            for refe in list_ref:
                shutil.move(refe, os.path.join(wdir, "library_cluster", refe.split('/')[-1]))

            shutil.copy(os.path.join(wdir, directory_database, 'list_rmsd.txt'),
                        os.path.join(wdir, "library_cluster", 'list_rmsd.txt'))
            tar_dir(os.path.join(wdir, directory_database), os.path.join(wdir, os.path.basename(directory_database)+'_lib.tar.gz'))

def get_artificial_parser_option_for_library_generation(local_grid=None,remote_grid=None,supercomputer=None,force_core=-1, directory_database=None,
                     c_angle=95, c_dist=95, c_angle_dist=95, c_cvl_diff=95, score_intra_fragment=95, j_angle=90, j_dist=90, j_angle_dist=90,
                     j_cvl_diff=95, score_inter_fragments=90, rmsd_clustering=1.5, exclude_residues_superpose=0, work_directory="./",
                     targz=None, pdbmodel=None, remove_coil=False, weight = "distance_avg",
                     strictness_ah=0.45, strictness_bs=0.45, peptide_length=3, enhance_fold=False, superpose=True, process_join=False,
                     nilges=3, sequence="", ncssearch=False, multimer=True,
                     rmsd_min=0.0, rmsd_max=6.0, ssbridge=False, connectivity=False, representative=False, sidechains=False,
                     doCluster=True, sym=None, criterium_selection_core="residues",
                     test=False, cath_id=None, target_sequence=None, swap_superposition=False,
                     clustering_mode='no_clustering', number_of_ranges=500, number_of_clusters=7000, exclude_sequence=None, reverse=False,use_model_as_it_is=False, superposition_mode=False):

    def inner(*args, **kwargs):
        pars = SystemUtility.AttrDict()
        pars.update(kwargs)
        return pars

    return inner(local_grid = local_grid, remote_grid = remote_grid, supercomputer = supercomputer, force_core = force_core, directory_database = directory_database,
    c_angle = c_angle, c_dist = c_dist, c_angle_dist = c_angle_dist, c_cvl_diff = c_cvl_diff, score_intra_fragment = score_intra_fragment, j_angle = j_angle, j_dist = j_dist, j_angle_dist = j_angle_dist,
    j_cvl_diff = j_cvl_diff, score_inter_fragments = score_inter_fragments, rmsd_clustering = rmsd_clustering, exclude_residues_superpose = exclude_residues_superpose, work_directory = work_directory,
    targz = targz, pdbmodel = pdbmodel, remove_coil = remove_coil, weight = weight,
    strictness_ah = strictness_ah, strictness_bs = strictness_bs, peptide_length = peptide_length, enhance_fold = enhance_fold, superpose = superpose, process_join = process_join,
    nilges = nilges, sequence=sequence, ncssearch=ncssearch, multimer=multimer,
    rmsd_min=rmsd_min, rmsd_max=rmsd_max, ssbridge=ssbridge, connectivity=connectivity, representative=representative, sidechains=sidechains,
    criterium_selection_core=criterium_selection_core, test=test, cath_id=cath_id, target_sequence=target_sequence,
    clustering_mode=clustering_mode, number_of_ranges=number_of_ranges, number_of_clusters=number_of_clusters, exclude_sequence=exclude_sequence,use_model_as_it_is=use_model_as_it_is, superposition_mode=superposition_mode)

@SystemUtility.timing
def generate_library(local_grid=None,remote_grid=None,supercomputer=None,force_core=-1, directory_database=None,
                     c_angle=95, c_dist=95, c_angle_dist=95, c_cvl_diff=95, score_intra_fragment=95, j_angle=90, j_dist=90, j_angle_dist=90,
                     j_cvl_diff=95, score_inter_fragments=90, rmsd_clustering=1.5, exclude_residues_superpose=0, work_directory="./",
                     targz=None, pdbmodel=None, remove_coil=False, weight = "distance_avg",
                     strictness_ah=0.45, strictness_bs=0.45, peptide_length=3, enhance_fold=False, superpose=True, process_join=False,
                     nilges=3, sequence="", ncssearch=False, multimer=True,
                     rmsd_min=0.0, rmsd_max=6.0, ssbridge=False, connectivity=False, representative=False, sidechains=False,
                     doCluster=True, sym=None, criterium_selection_core="residues",
                     test=False, cath_id=None, target_sequence=None, swap_superposition=False,
                     clustering_mode='no_clustering', number_of_ranges=500, number_of_clusters=7000, exclude_sequence=None, reverse=False,
                     use_model_as_it_is=False, superposition_mode=False, interface=False, verbose=False):

    global PATH_NEW_BORGES
    global GRID_TYPE_R
    global MAX_PDB_TAR
    global toclientdir
    global doCluster_global
    global PYTHON_V
    global PATH_PYTHON_INTERPRETER

    MAX_NUM_FOR_TEST = 100
    MIN_NON_RANDOM_SCORE = 45.0
    if swap_superposition:
        criterium_selection_core += "|||swap"
        if reverse:
            criterium_selection_core += "reverse"

    if sym is None:
        sym = SystemUtility.SystemUtility()

    if local_grid is not None or remote_grid is not None or supercomputer is not None:
        SystemUtility.startCheckQueue(sym, delete_check_file=False, callback=evaluate_grid_job, forcework=True)

    if score_intra_fragment < MIN_NON_RANDOM_SCORE: score_intra_fragment = MIN_NON_RANDOM_SCORE
    if score_inter_fragments < MIN_NON_RANDOM_SCORE: score_inter_fragments = MIN_NON_RANDOM_SCORE

    if int(force_core) > 0:
        sym.PROCESSES = int(force_core)
    if c_angle <= 0: #and directory_database is not None:
        c_angle = score_intra_fragment
    if c_dist <= 0: #and directory_database is not None:
        c_dist = score_intra_fragment
    if c_angle_dist <= 0: # and directory_database is not None:
        c_angle_dist = score_intra_fragment
    if c_cvl_diff <= 0: # and directory_database is not None:
        c_cvl_diff = score_intra_fragment
    if j_angle <= 0: # and directory_database is not None:
        j_angle = score_inter_fragments
    if j_dist <= 0: # and directory_database is not None:
        j_dist = score_inter_fragments
    if j_angle_dist <= 0: # and directory_database is not None:
        j_angle_dist = score_inter_fragments
    if j_cvl_diff <= 0: # and directory_database is not None:
        j_cvl_diff = score_inter_fragments

    if interface:
        multimer = False
        # representative = True TODO: maybe it is better to active when you are using representative (?)

    thresh = float(rmsd_clustering)
    if thresh <= 0.0 or clustering_mode == 'no_clustering' or clustering_mode == 'rmsd_range' or clustering_mode == 'random_sampling':
        doCluster = False

    doCluster_global = doCluster

    superpose_exclude = int(exclude_residues_superpose) + 1
    strucc = None
    pattern_cvs = None
    pattern = None
    lisBig = None
    highd = None
    wdir = work_directory
    cm = None
    toclientdir = None
    pdbsearch = None

    if wdir is None or wdir == "":
        wdir = "./"

    if not os.path.exists(wdir):
        os.makedirs(wdir)


    DicGridConn = {}
    baseline = []
    skip = False
    bl = sys.argv[1:]
    for arg in bl:
        if skip:
            skip = False
            continue
        if arg in ["--directory_database", "--supercomputer", "--local_grid", "--remote_grid", "--pdbmodel"]:
            skip = True
            continue
        baseline.append(arg)

    if target_sequence is not None:
        tg = ""
        with open(target_sequence, "r") as f:
            for line in f.readlines():
                if not line.startswith(">"): tg += line[:-1]

        hits = get_blast_models(tg)
        top5 = [(hits[key]['pdbid'],hits[key]['chains']) for key in sorted(hits.keys(), key=lambda x: hits[x]["evalue"])][:5]
        pfams,cathids,scops = get_structure_ids_from_pdbids(top5)
        directory_database = download_database(pfams=pfams,cathids=cathids,scops=scops)
        if directory_database is None:
            print("ERROR: Impossible to generate the database from sequence")
            raise Exception("ERROR: Impossible to generate the database from sequence")
    elif cath_id is not None:
        directory_database = download_database(cathids=[cath_id])
        if directory_database is None:
            print("ERROR: Impossible to generate the database from the CATHID")
            raise Exception("ERROR: Impossible to generate the database from the CATHID")

    if (local_grid is not None and os.path.exists(local_grid)) or (remote_grid is not None and os.path.exists(remote_grid)):
        if targz is not None:
            print(colored("ATTENTION: --targz option is incompatible with --remote_grid , --local_grid and --supercomputer", "red"))
            print()
            #traceback.print_exc(file=sys.stdout)
            raise Exception("--targz option is incompatible with --remote_grid , --local_grid and --supercomputer")

        toclientdir = os.path.join(wdir, "./ToProcess")

        if not os.path.exists(toclientdir):
            os.makedirs(toclientdir)

        # read the setup.bor and configure the grid
        setupbor = None
        GRID_TYPE = None
        if remote_grid is not None and os.path.exists(remote_grid):
            path_bor = remote_grid
            try:
                if PYTHON_V == 3:
                    setupbor = configparser.ConfigParser()
                    setupbor.read_file(open(path_bor))
                elif PYTHON_V == 2:
                    setupbor = ConfigParser.ConfigParser()
                    setupbor.readfp(open(path_bor))
                DicGridConn["username"] = setupbor.get("GRID", "remote_frontend_username")
                DicGridConn["host"] = setupbor.get("GRID", "remote_frontend_host")
                DicGridConn["port"] = setupbor.getint("GRID", "remote_frontend_port")
                DicGridConn["passkey"] = setupbor.get("CONNECTION", "remote_frontend_passkey")
                DicGridConn["promptA"] = (setupbor.get("GRID", "remote_frontend_prompt")).strip() + " "
                DicGridConn["isnfs"] = setupbor.getboolean("GRID", "remote_fylesystem_isnfs")
                try:
                    DicGridConn["remote_submitter_username"] = setupbor.get("GRID", "remote_submitter_username")
                    DicGridConn["remote_submitter_host"] = setupbor.get("GRID", "remote_submitter_host")
                    DicGridConn["remote_submitter_port"] = setupbor.getint("GRID", "remote_submitter_port")
                    DicGridConn["promptB"] = (setupbor.get("GRID", "remote_submitter_prompt")).strip() + " "
                except:
                    pass
                DicGridConn["home_frontend_directory"] = setupbor.get("GRID", "home_frontend_directory")
                PATH_NEW_BORGES = setupbor.get("GRID", "path_remote_borges")
                PATH_PYTHON_INTERPRETER = setupbor.get("GRID", "python_remote_interpreter")
                if PATH_PYTHON_INTERPRETER.strip() in ["", None]: PATH_PYTHON_INTERPRETER = "/usr/bin/python"
                GRID_TYPE = setupbor.get("GRID", "type_remote")
            except:
                print(colored("ATTENTION: Some keyword in your configuration files are missing. Contact your administrator", "red"))
                print("Path bor given: ", path_bor)
                print()
                traceback.print_exc(file=sys.stdout)
                raise Exception("Some keyword in your configuration files are missing. Contact your administrator")
        elif local_grid is not None and os.path.exists(local_grid):
            path_bor = local_grid
            try:
                if PYTHON_V == 3:
                    setupbor = configparser.ConfigParser()
                    setupbor.read_file(open(path_bor))
                elif PYTHON_V == 2:
                    setupbor = ConfigParser.ConfigParser()
                PATH_NEW_BORGES = setupbor.get("LOCAL", "path_local_borges")
                PATH_PYTHON_INTERPRETER = setupbor.get("LOCAL", "python_local_interpreter")
                if PATH_PYTHON_INTERPRETER.strip() in ["", None]: PATH_PYTHON_INTERPRETER = "/usr/bin/python"
                GRID_TYPE = setupbor.get("GRID", "type_local")
            except:
                print(colored("ATTENTION: Some keyword in your configuration files are missing. Contact your administrator","red"))
                print("Path bor given: ", path_bor)
                print()
                traceback.print_exc(file=sys.stdout)
                raise Exception("Some keyword in your configuration files are missing. Contact your administrator")

        # STARTING THE GRID MANAGER
        if cm is None:
            if GRID_TYPE == "Condor":
                cm = Grid.condorManager()
            elif GRID_TYPE == "SGE":
                QNAME = setupbor.get("SGE", "qname")
                FRACTION = setupbor.getfloat("SGE", "fraction")
                cm = Grid.SGEManager(qname=QNAME, fraction=FRACTION)
            elif GRID_TYPE == "MOAB":
                PARTITION = setupbor.get("MOAB", "partition")
                # FRACTION = setupbor.getfloat("MOAB","partition")
                cm = Grid.MOABManager(partition=PARTITION)
            elif GRID_TYPE == "SLURM":
                PARTITION = setupbor.get("SLURM", "partition")
                if PARTITION is not None and PARTITION != '':
                    cm = Grid.SLURMManager(partition=PARTITION)
                else:
                    cm = Grid.SLURMManager()
            elif GRID_TYPE == "TORQUE":
                QNAME = setupbor.get("TORQUE", "qname")
                FRACTION = setupbor.getint("TORQUE", "cores_per_node")
                PARALLEL_JOBS = setupbor.getint("TORQUE", "number_of_parallel_jobs")
                MAUI = setupbor.getboolean("TORQUE", "maui")
                cm = Grid.TORQUEManager(qname=QNAME, cores_per_node=FRACTION, parallel_jobs=PARALLEL_JOBS, maui=MAUI)

        if cm is not None:
            cm.setRank("kflops")
            cm.nice_user = "true"

    elif supercomputer is not None and os.path.exists(supercomputer):
        # read the nodes.txt file and process it
        #NOTE:TODO: It must be read also created the PATH_PYTHON_INTERPRETER by extracting this info from nodelist which should be created accordingly
        toclientdir = os.path.join(wdir, "./ToProcess")

        if not os.path.exists(toclientdir):
            os.makedirs(toclientdir)

        path_nodes = supercomputer
        f = open(path_nodes, "r")
        nodes_list = f.readlines()
        f.close()
        PATH_NEW_BORGES = nodes_list[0]
        nodes_list = nodes_list[1:]
        for i in range(len(nodes_list)):
            nodes_list[i] = nodes_list[i][:-1] + "***" + str(i)
        SystemUtility.NODES = nodes_list

    if targz is None and directory_database is not None:
        pdbsearchin = ""
        listn = []
        if not use_model_as_it_is:
            graph_full_reference, stru, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(pdbmodel,
                                                                                                                         weight=weight,
                                                                                                                         strictness_ah=strictness_ah,
                                                                                                                         strictness_bs=strictness_bs,
                                                                                                                         peptide_length=peptide_length,
                                                                                                                         write_pdb=False, only_reformat=True)
            all_frags = get_all_fragments(graph_full_reference)

            if remove_coil:
                for fra in all_frags:
                    for resi in fra["resIdList"]:
                        listn.append((fra["chain"], resi))

            for model in stru.get_list():
                reference = []
                for chain in model.get_list():
                    for residue in chain.get_list():
                        if len(listn) == 0 or (chain.get_id(), residue.get_id()) in listn:
                            reference += residue.get_unpacked_list()

                pdbmod, cnv = Bioinformatics.get_pdb_from_list_of_atoms(reference, renumber=True, uniqueChain=False, remove_non_res_hetatm=True, polyala=not sidechains)
                # pdbmod = "MODEL "+str(model.get_id())+"\n"+pdbmod+"\n\n"
                pdbsearchin += pdbmod
                # pdbsearchin += "ENDMDL\n\n"

            fds = open(os.path.join(wdir, os.path.basename(pdbmodel)[:-4] + "_input_search.pdb"), "w")
            fds.write(pdbsearchin)
            fds.close()
        else:
            shutil.copyfile(pdbmodel,os.path.abspath(os.path.join(wdir, os.path.basename(pdbmodel)[:-4] + "_input_search.pdb")))
            with open(pdbmodel,"r") as f: pdbsearchin = f.read()

        #Forcing the garbage collector to empty space
        listn = []
        tupls = []

        pdbmodel = os.path.abspath(os.path.join(wdir, os.path.basename(pdbmodel)[:-4] + "_input_search.pdb"))

    if (local_grid is not None and os.path.exists(local_grid)) or (
                    remote_grid is not None and os.path.exists(remote_grid)) or (
                    supercomputer is not None and os.path.exists(supercomputer)):
        shutil.copyfile(pdbmodel, os.path.join(toclientdir, os.path.basename(pdbmodel)))
        baseline.append("--pdbmodel")
        baseline.append(os.path.basename(pdbmodel))

    pattern = None
    pattern_cvs = None
    highd = None
    if targz is not None and os.path.exists(targz):
        fileParaName = targz
        t = datetime.datetime.now()
        epcSec = time.mktime(t.timetuple())
        now = datetime.datetime.fromtimestamp(epcSec)
        print("" + str(now.ctime()))
        sys.stdout.flush()

        tar = tarfile.open(fileParaName, "r:gz")
        print("Read 1 Parameter: the tar.gz archive")
        infile = tar.extractfile((fileParaName[:-7]) + "PARAM")
        print("Read 2 Parameter: the number of pdbs to work with")
        nPDBs = pickle.load(infile)
        pdbsol = []
        solred = {}

        dic_json = {'number of targets': len(nPDBs)}
        modify_json_file('output.json', dic_json, lib_generation=True)

        for i in range(nPDBs):
            print("Read 3 Parameter: PDB file to work with")
            pdb_name = pickle.load(infile)
            pdbfile = pickle.load(infile)
            cvs_list_str = pickle.load(infile)
            pdbstruc = pickle.load(infile)
            print ("pdbf",pdb_name)
            if i == 0:
                pattern, pattern_cvs, highd, pdbsearch, strucc, lisBig, peptide_length = evaluate_model(pdbmodel, bool(enhance_fold), int(peptide_length),  weight, strictness_ah, strictness_bs)
            # pdbsearch = pdbsearchin
            pdl, dicl = evaluate_pdb(sym, wdir, pdbfile, cvs_list_str, pdbstruc, strucc, lisBig, i, pattern,
                                     pattern_cvs, highd, doCluster, superpose, process_join, pdbsearch, pdb_name,
                                     thresh, superpose_exclude, peptide_length, sequence,
                                     ncssearch, multimer, c_angle, c_dist, c_angle_dist, c_cvl_diff, j_angle, j_dist,
                                     j_angle_dist, j_cvl_diff, rmsd_min, rmsd_max, ssbridge, connectivity,
                                     nilges, enhance_fold, representative, pdbmodel, sidechains, weight, strictness_ah,
                                     strictness_bs, criterium_selection_core, exclude_sequence, superposition_mode,
                                     return_pdbstring=True, interface=interface, verbose=verbose)
            pdbsol += pdl
            solred.update(dicl)
        infile.close()
        tar.close()

        fileout = os.path.basename(fileParaName[:-7]) + "_res"
        pdbf = os.path.join(wdir, fileout + ".tar")
        fro = open(os.path.join(wdir, fileout + "_out.data"), "wb")
        # print "Write number of solutions...",len(pdbsol)
        # pickle.dump(len(pdbsol),fro)
        pickle.dump(float(thresh), fro)
        pickle.dump(int(superpose_exclude), fro)
        pickle.dump(int(nilges), fro)
        pickle.dump(weight, fro)
        pickle.dump(strictness_ah, fro)
        pickle.dump(strictness_bs, fro)
        pickle.dump(peptide_length, fro)
        pickle.dump(criterium_selection_core, fro)
        pickle.dump(pdbsol, fro)
        pickle.dump(solred, fro)
        fro.close()
        tar = tarfile.open(pdbf, "a")
        tar.add(os.path.join(wdir, fileout + "_out.data"), arcname=fileout + "_out.data")
        tar.close()
        os.remove(os.path.join(wdir, fileout + "_out.data"))
        compri = gzip.open(pdbf + ".gz", 'wb')
        fion = open(pdbf, "rb")
        compri.write(fion.read())
        fion.close()
        compri.close()
        os.remove(pdbf)
    elif directory_database is not None:
        allfiles = []
        if os.path.exists(os.path.join(directory_database, "listfiles.txt")):
            print("Reading", os.path.join(directory_database, "listfiles.txt"), "...")
            f = open(os.path.join(directory_database, "listfiles.txt"), "r")
            alllinesaa = f.readlines()
            f.close()
            print("Done...")
            for pdbfl in alllinesaa:
                if len(allfiles) % 1000 == 0:
                    print("Parsed", len(allfiles))
                pdbfl = pdbfl.split()
                pdbf = pdbfl[0]
                if pdbf.endswith(".pdb") or pdbf.endswith(".gz") or pdbf.endswith(".ent"):
                    allfiles.append(pdbf)
                elif pdbf.endswith(".data"):
                    for key in pdbfl[1:]:
                        allfiles.append((pdbf, key))
        else:
            for root, subFolders, files in os.walk(directory_database):
                for fileu in files:
                    pdbf = os.path.join(root, fileu)
                    if pdbf.endswith(".pdb") or pdbf.endswith(".gz") or pdbf.endswith(".ent"):
                        allfiles.append(pdbf)
                    elif pdbf.endswith(".data"):
                        cv_matrices_shelve = shelve.open(pdbf)
                        for key in cv_matrices_shelve:
                            allfiles.append((pdbf, key))

        pdbsfiles = []
        print("Starting shuffle...")
        random.shuffle(allfiles)
        print("List shuffled!")
        if test: allfiles = allfiles[:MAX_NUM_FOR_TEST]

        dic_json = {'number of targets': len(allfiles)}
        modify_json_file('output.json', dic_json, lib_generation=True)

        if interface:
            dir_spack = os.path.join(os.getcwd(), 'dir_symmetry')
            if os.path.exists(dir_spack):
                shutil.rmtree(dir_spack)
            os.mkdir(dir_spack)

        for i in range(len(allfiles)):
            pdbf = allfiles[i]
            if isinstance(pdbf, str):
                if pdbf.endswith(".gz"):
                    fileObj = gzip.GzipFile(pdbf, 'rb')
                    fileContent = fileObj.read()
                    fileObj.close()
                    os.remove(pdbf)
                    pdbf = pdbf[:-3]  # elimino estensione .gz
                    fou = open(pdbf, "w")
                    fou.write(fileContent.decode("utf-8") if isinstance(fileContent,bytes) else fileContent)
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
                    Bioinformatics.rename_hetatm_and_icode(pdbf)

                if interface:
                    pdbf_copy = os.path.join(os.getcwd(), os.path.basename(pdbf))
                    shutil.copyfile(pdbf, pdbf_copy)
                    crystcard = Bioinformatics.get_crystcard(pdbf_copy)
                    input_spack = Bioinformatics.prepare_model_for_spack(pdbf_copy, crystcard)
                    distance = 10  # TODO Change it in the future
                    generate_symmetry_copies(os.path.basename(input_spack), radius=str(distance))
                    output_spack = input_spack[:-4] + 'sym.pdb'
                    output_spack_new = os.path.join(dir_spack, os.path.basename(output_spack)[:-4] + str(distance) + '.pdb')
                    os.rename(output_spack, output_spack_new)
                    pdbf = output_spack_new
                    os.remove(input_spack)
                    os.remove(pdbf_copy)

                    # crystcard = Bioinformatics.get_crystcard(pdbf)
                    # graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
                    #     pdbf, weight=weight, strictness_ah=strictness_ah, strictness_bs=strictness_bs,
                    #     peptide_length=peptide_length,write_pdb=True)
                    # model_annotated = os.path.join(os.path.basename(pdbf[:-4] + '_input_search.pdb'))
                    # input_spack = Bioinformatics.prepare_model_for_spack(model_annotated, crystcard)
                    # distance = 10 #TODO Change it in the future
                    # generate_symmetry_copies(input_spack, radius=str(distance))
                    # output_spack = input_spack[:-4] + 'sym.pdb'
                    # # output_spack_new = os.path.join(dir_spack, os.path.basename(pdbf)[:-4] + '_sym' + str(distance) + '.pdb')
                    # output_spack_new = os.path.join(dir_spack, output_spack[:-4] + str(distance) + '.pdb')
                    # os.rename('./' + output_spack, output_spack_new)
                    # os.remove(model_annotated)
                    # os.remove(input_spack)
                    # pdbf = output_spack_new
            if i == 0:
                pattern, pattern_cvs, highd, pdbsearch, strucc, lisBig, peptide_length = evaluate_model(pdbmodel, bool(enhance_fold), int(peptide_length),  weight, strictness_ah, strictness_bs)
                pdbsearch = pdbsearchin
                try:
                    shutil.copyfile(pdbmodel,os.path.join(toclientdir, os.path.basename(pdbmodel)))
                except:
                    pass

            # Multiprocessing
            if supercomputer is None and local_grid is None and remote_grid is None:
                if isinstance(pdbf, str):
                    evaluate_pdb(sym, wdir, pdbf, None, None, strucc, lisBig, i, pattern, pattern_cvs,
                                 highd, doCluster, superpose, process_join, pdbsearch, "", thresh, superpose_exclude,
                                 peptide_length,  sequence, ncssearch, multimer, c_angle, c_dist, c_angle_dist, c_cvl_diff, j_angle, j_dist,
                                     j_angle_dist, j_cvl_diff, rmsd_min, rmsd_max, ssbridge, connectivity,
                                     nilges, enhance_fold, representative, pdbmodel, sidechains, weight, strictness_ah,
                                 strictness_bs, criterium_selection_core,
                                 exclude_sequence, superposition_mode, interface=interface, verbose=verbose)

                elif isinstance(pdbf, tuple):
                    pdbpa = pdbf[0]
                    key = pdbf[1]
                    cv_matrices_shelve = shelve.open(pdbpa)
                    cvs_list_str = cv_matrices_shelve[key]["cvs_list"]
                    pdbstruc = cv_matrices_shelve[key]["structure"]
                    evaluate_pdb(sym, wdir, pdbpa, cvs_list_str, pdbstruc, strucc, lisBig, i, pattern,
                                 pattern_cvs, highd, doCluster, superpose, process_join, pdbsearch, key, thresh,
                                 superpose_exclude, peptide_length, sequence,
                                     ncssearch, multimer, c_angle, c_dist, c_angle_dist, c_cvl_diff, j_angle, j_dist,
                                     j_angle_dist, j_cvl_diff, rmsd_min, rmsd_max, ssbridge, connectivity,
                                     nilges, enhance_fold, representative, pdbmodel, sidechains, weight,
                                 strictness_ah, strictness_bs, criterium_selection_core,
                                 exclude_sequence, superposition_mode, interface=interface, verbose=verbose)

            else:
                if isinstance(pdbf, str):
                    f = open(pdbf, "r")
                    pdbread = f.read()
                    f.close()
                    pdbsfiles.append((os.path.basename(pdbf)[:-4], pdbread, None, None))
                elif isinstance(pdbf, tuple):
                    pdbpa = pdbf[0]
                    key = pdbf[1]
                    cv_matrices_shelve = shelve.open(pdbpa)
                    cvs_list_str = cv_matrices_shelve[key]["cvs_list"]
                    pdbstruc = cv_matrices_shelve[key]["structure"]
                    pdbsfiles.append((key, key, cvs_list_str, pdbstruc))

                if (i != 0 and i % MAX_PDB_TAR == 0) or i == len(allfiles) - 1:
                    prepare_and_launch_job(cm, baseline, "job_" + str(i), pdbsfiles, supercomputer, pdbmodel)
                    pdbsfiles = []
            if perform_superposition and superposition_mode:
                print("Best superposition saved as", os.path.basename(pdbmodel)[:-17] + '.pdb')

    if local_grid is not None or remote_grid is not None or supercomputer is not None:
        SystemUtility.endCheckQueue()
    print("...Writing output files and cleaning...")

    cluster_library(directory_database="./library", rmsd_clustering=rmsd_clustering, clustering_mode=clustering_mode, number_of_ranges=number_of_ranges,
                    number_of_clusters=number_of_clusters, exclude_residues_superpose=exclude_residues_superpose, ssbridge=ssbridge, nilges=nilges,
                    work_directory=wdir)

def generate_symmetry_copies (model, radius='10'):

    platform = sys.platform
    if platform == 'darwin':
        spack = os.path.join(os.path.dirname(sys.argv[0]), 'Spack_mac') #Relative path to spack for mac
    elif platform == 'linux' or platform == 'linux2':
        spack = os.path.join(os.path.dirname(sys.argv[0]), 'Spack_linux')  # Relative path to spack for linux
    else:
        print("ERROR: The operating system is not compatible with the Spack version")
        raise Exception("The operating system is not compatible with the Spack version")

    if not os.path.isfile(spack):
        print("ERROR: Spack software cannot be found")
        raise Exception("Spack software cannot be found")

    dist = '-D' + radius
    args_spack = [spack, '-A', dist, model]

    # f = open('output_spack.log', 'w') # Uncomment if you want to capture log
    # p = subprocess.Popen(args_spack, stdout=f, stderr=f) # Uncomment if you want to capture log
    # f.close()  # Uncomment if you want to capture log
    p = subprocess.Popen(args_spack,  shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE) # Comment if you want to capture log
    out, err = p.communicate()  # Finishing subprocess

def write_from_graph(graph, strucc, path_out, name_out, write_pdb=True, write_graphml=False, interfaces = False):

    g = copy.deepcopy(graph) #For preventing overwritting although it is not necessary

    if write_pdb:
        pdbsearchin = ''
        listar = []

        for frag in g.vs:
            for res in frag["reslist"]:
                if len(res) == 5:
                    listar.append(tuple(res[:-1]))
                else:
                    listar.append(tuple(res))

        for model in strucc.get_list():
            reference = [atm for atm in Bio.PDB.Selection.unfold_entities(strucc[model.get_id()], "A") if (len(listar) == 0 or atm.get_parent().get_full_id() in listar)]

            if interfaces:
                pdbString, cnv = Bioinformatics.get_pdb_from_list_of_atoms(reference, renumber=False, uniqueChain=False,remove_non_res_hetatm=False, sort_reference=False, chainFragment=True, interface=True)
                pdbsearchin += pdbString

            else:
                try:
                    pdbsearchin = Bioinformatics.generate_secondary_structure_record(g)
                except:
                    print("No se puede generar el HELIX/SHEET record del pdb")
                pdbString, cnv = Bioinformatics.get_pdb_from_list_of_atoms(reference, renumber=False, uniqueChain=False,
                                                                           remove_non_res_hetatm=False,
                                                                           sort_reference=False, chainFragment=False,
                                                                           interface=True) #NOTE AMB: I have change chainFragment to False

                pdbsearchin += pdbString
        fds = open(os.path.join(path_out, os.path.basename(name_out)) + ".pdb", "w")
        fds.write(pdbsearchin)
        fds.close()
        print("The pdb file has been writen:", os.path.join(path_out, os.path.basename(name_out)) + ".pdb")
    if write_graphml:
        g.write_graphml(os.path.join(path_out, os.path.basename(name_out)) + ".graphml")
        print("The graphml file has been writen:", os.path.join(path_out, os.path.basename(name_out)) + ".graphml")

def generate_graph_interfaces (dict_pdbs, get_unique=True):

    graph_interfaces = Graph()
    number_node = 0

    for chains in dict_pdbs:
        last_node = None
        for ind, distance in enumerate(dict_pdbs[chains]):
            if get_unique == True and dict_pdbs[chains][distance]['unique'] == 0: #Writing just unique interfaces
                continue

            graph_interfaces.add_vertices(1)
            graph_interfaces.vs[number_node]["pdb"] = dict_pdbs[chains][distance]['pdbf']
            graph_interfaces.vs[number_node]["interface"] = chains
            graph_interfaces.vs[number_node]["avg_distance"] = distance
            # graph_interfaces.vs[number_node]["unique"] = dict_pdbs[chains][distance]['unique']
            graph_interfaces.vs[number_node]["number_frags"] = len(dict_pdbs[chains][distance]['graph'].vs())
            graph_interfaces.vs[number_node]["reslist"] = [ind['reslist'] for ind in dict_pdbs[chains][distance]['graph'].vs()]
            graph_interfaces.vs[number_node]["sequence"] = [ind['sequence'] for ind in dict_pdbs[chains][distance]['graph'].vs()]
            graph_interfaces.vs[number_node]["number_res"] = [len(ind['sequence']) for ind in dict_pdbs[chains][distance]['graph'].vs()]
            graph_interfaces.vs[number_node]["avg_bfact"] = [ind['avg_bfact'] for ind in dict_pdbs[chains][distance]['graph'].vs()]
            graph_interfaces.vs[number_node]["centrality"] = [ind['centrality'] for ind in dict_pdbs[chains][distance]['graph'].vs()]

            if last_node != None:
                graph_interfaces.add_edges([(number_node, last_node)])
                graph_interfaces.es[graph_interfaces.get_eid(number_node, last_node)]["min"] = 17 #TODO: propiedad que añadir al edge


            last_node = number_node
            number_node += 1

    return(graph_interfaces)

@SystemUtility.timing
def analyse_interfaces(reference, strictness_ah, strictness_bs, peptide_length, max_ah_dist, min_ah_dist, max_bs_dist,
                       min_bs_dist, write_graphml, write_pdb, min_distance=10, max_distance=20, edges_distance=20,
                       verbose=False):

    if strictness_ah < 0:
        strictness_ah = 0.45
    if strictness_bs < 0:
        strictness_bs = 0.20
    if min_bs_dist < 0:
        min_bs_dist = 0
    if max_bs_dist < 0:
        max_bs_dist = 15
    if min_ah_dist < 0:
        min_ah_dist = 0
    if max_ah_dist < 0:
        max_ah_dist = 20
    if min_distance and min_distance < 1:
        min_distance = 10
    if max_distance and max_distance < min_distance:
        max_distance = min_distance + 1
    if edges_distance < 0:
        edges_distance = 20

    weight = "distance_avg"
    pep_len = int(peptide_length)

    crystcard = Bioinformatics.get_crystcard(reference)
    graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
        reference, weight=weight, strictness_ah=strictness_ah, strictness_bs=strictness_bs, peptide_length=pep_len,
        write_pdb=write_pdb, interface = False)
    graph_no_coil_reference = graph_full_reference.vs.select(sstype_in=["ah", "bs"]).subgraph()
    for frag in graph_no_coil_reference.vs(): frag['chain'] = frag['reslist'][0][2] #Adding chain to vertex/node
    for edge in graph_no_coil_reference.es(): edge['chains'] = [graph_no_coil_reference.vs[edge.source].attributes()['chain'],
                                                                graph_no_coil_reference.vs[edge.target].attributes()['chain']] #Adding chain to edges

    # write_from_graph(graph_no_coil_reference, structure_reference, './', reference[:-4] + '_jaja', write_pdb=False, write_graphml=True)
    list_distance = [edge['avg'] for edge in graph_no_coil_reference.es() if edge['chains'][0] == edge['chains'][1]] #Selecting distances withing the monomer
    # list_distance = [edge['avg'] for edge in graph_no_coil_reference.es()] #In case we want the max distance from the complete structure
    max_dist = max(list_distance)
    print("Max distances between fragments within the same monomer", max_dist)
    print("25% of max distance is", max_dist*0.25)
    small_radius = round((max_dist*0.25))
    big_radius = round((max_dist*0.4))
    print("Small radius is", small_radius)
    print("Big radius is", big_radius)


    if min_distance:
        small_radius = min_distance
    if max_distance:
        big_radius = max_distance

    model_annotated = os.path.basename(reference)
    # model_annotated = os.path.join(os.path.basename(reference[:-4] + '_input_search.pdb')) #In case the annotated model is used for spack
    input_spack = Bioinformatics.prepare_model_for_spack (model_annotated, crystcard)

    dir_spack = os.path.join(os.getcwd(), 'dir_symmetry')
    if not os.path.exists(dir_spack):
        os.mkdir(dir_spack)

    generate_symmetry_copies(input_spack, radius=str(small_radius))
    output_spack = input_spack[:-4] + 'sym.pdb'
    output_spack_new = os.path.join(dir_spack, output_spack[:-4] + str(small_radius) + '.pdb')
    os.rename('./' + output_spack, output_spack_new)
    os.remove(input_spack)

    #GENERATING THE GRAPH
    graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
            output_spack_new, weight=weight, strictness_ah=strictness_ah, strictness_bs=strictness_bs, peptide_length=pep_len,
            write_pdb=False, min_ah=1, min_bs=1, interface = True, verbose=False)
    graph_no_coil_reference = graph_full_reference.vs.select(sstype_in=["ah", "bs"]).subgraph() #Removing coils until spack's bug is gone
    write_from_graph(graph_no_coil_reference, structure_reference, './', output_spack[:-4] + str(small_radius) + '_input_search', write_pdb=True, write_graphml=False)

    #ADDING AVG-"CA"-BFACTORS AND CHAIN-ID TO EACH NODE
    for i, frag in enumerate(graph_no_coil_reference.vs()):
        frag['chain'] = frag['reslist'][0][2] + frag['reslist'][0][3][2] # Chain = chain + insertion_code (represents the sym operation)
        frag['avg_bfact'] = 0
        for residue in frag['reslist']:
            model = residue[1]
            chain = residue[2][0]
            residue = residue[3]
            atom = structure_reference[model][chain][residue]["CA"]
            frag['avg_bfact'] += atom.get_bfactor()
        frag['avg_bfact'] = round(frag['avg_bfact']/len(frag['reslist']), 3)

    #ADDING CHAIN ID FROM SOURCE AND TARGET VERTEX TO EDGES
    for edge in graph_no_coil_reference.es():
        chain_source = graph_no_coil_reference.vs[edge.source].attributes()['chain']
        chain_target = graph_no_coil_reference.vs[edge.target].attributes()['chain']
        edge['chains'] = [chain_source, chain_target]


    # COMPUTING CENTRALITY #This is not used right now
    edges_filter = graph_no_coil_reference.es.select(min_le=small_radius)
    graph_no_coil_filter = graph_no_coil_reference.subgraph_edges(edges_filter,delete_vertices=False)

    eigen_ref = graph_no_coil_filter.evcent(directed=False, scale=True,
                                            weights=graph_no_coil_filter.es["weight"],
                                            return_eigenvalue=True)
    # print("EIGEN Centrality values:", eigen_ref, len(eigen_ref[0]), len(graph_no_coil_filter.vs))

    for i, frag in enumerate(graph_no_coil_reference.vs):
        frag['centrality'] = eigen_ref[0][i]



    #GETTING THE COMBINATIONS OF 2 CHAINS ORDERED BY NUMBER OF FRAGMENTS
    dict_chains = {}

    for node in graph_no_coil_reference.vs():
        if node['chain'] not in dict_chains:
            dict_chains[node['chain']] = 1
        else:
            dict_chains[node['chain']] += 1

    ordered_chains = [k for k, v in sorted(dict_chains.items(), key=lambda item: item[1], reverse=True)]

    # l = list(permutations(ordered_chains, 2)) #Permutations
    l = list(combinations(ordered_chains, 2)) #Combinations
    print("Combinations of chains ordered by number of fragments in them is:", l)

    dir_interfaces = os.path.join(os.getcwd(), 'interfaces')
    if os.path.exists(dir_interfaces):
        shutil.rmtree(dir_interfaces)
    os.mkdir(dir_interfaces)
    # write_from_graph(graph_no_coil_reference, structure_reference, dir_interfaces, reference[:-4] + '_interface_full', write_pdb=True, write_graphml=True)

    dict_pdbs = {}
    to_delete = []

    for ind, perm in enumerate(l):
        graph_perm = graph_no_coil_reference.es.select(chains_in=[list(perm), list(perm)[::-1]]).subgraph()
        # write_from_graph(graph_perm, structure_reference, dir_interfaces,reference[:-4] + '_interface_' + str(perm[0]) + str(perm[1]), write_pdb=True, write_graphml=True)

        list_distances = list(range(int(small_radius),int(big_radius)))
        print("List_distances is", list_distances)

        for ind2, dist in enumerate(list_distances):
            print("\nEvaluating distance", dist)
            list_dist = [edge_perm['avg'] for edge_perm in graph_perm.es()]
            if min(list_dist) > dist:
                print('The distance between chains', str(perm[0]) + '-' + str(perm[1]), 'is higher than the minimum distance', str(dist)+'Å')
                continue

            graph_perm_filter = graph_perm.es.select(avg_le=dist).subgraph() #Taking the avg distance
            # graph_perm_filter = graph_perm.es.select(min_le=dist).subgraph() #Taking the min distance

            write_from_graph(graph_perm_filter, structure_reference, dir_interfaces,
                         reference[:-4] + '_interface_' + str(perm[0]) + str(perm[1]) + '_' + str(dist), write_pdb=True, write_graphml=False, interfaces = True)

            if not perm in dict_pdbs:
                dict_pdbs[perm] = {}
                dict_pdbs[perm][dist] = {}
                print(dict_pdbs[perm])
                dict_pdbs[perm][dist]['unique'] = 1

            dict_pdbs[perm][dist] = {}
            dict_pdbs[perm][dist]['graph'] = graph_perm_filter

            ########################## PDBFILE ##########################
            dict_pdbs[perm][dist]['pdbf'] = os.path.join(dir_interfaces,
                                                             os.path.basename(reference)[:-4] + '_interface_' + str(
                                                                 perm[0]) + str(perm[1]) + '_' + str(dist) + '.pdb')

            ########################## PDBSEARCH ##########################
            with open(dict_pdbs[perm][dist]['pdbf'], "r") as f: dict_pdbs[perm][dist]['pdbsearch'] = f.read()

            ########################## CVS_MODEL ##########################
            strucc, cvs_list = parse_pdb_and_generate_cvs(dict_pdbs[perm][dist]['pdbf'], strucc=None,
                                                              peptide_length=3, one_model_per_nmr=True,
                                                              only_reformat=True, process_only_chains=None)
            dict_pdbs[perm][dist]['cvs_model'] = [cv for frag in cvs_list for cv in frag]

            ########################## IDENTITY ##########################
            identity = (95.0, 95.0, 95.0, 95.0, 90.0, 90.0, 90.0, 90.0)

            list_chains = list(perm)

            dict_pdbs[perm][dist]['unique'] = 3  # Tag for indetermined interfaces

            if ind2 > 0 and list_distances[ind2-1] in dict_pdbs[perm]:
                print("Comparing", dict_pdbs[perm][list_distances[ind2]]['pdbf'])
                print("with", dict_pdbs[perm][list_distances[ind2-1]]['pdbf'])
                structure, solutions, chains_dict, lisBigSS = process_structure_fast(pdbsearch = dict_pdbs[perm][dist]['pdbsearch'],
                                                                                 pdbf= dict_pdbs[perm][list_distances[ind2-1]]['pdbf'],
                                                                                 pattern= None,
                                                                                 cvs_model = dict_pdbs[perm][dist]['cvs_model'],
                                                                                 sequence = '',
                                                                                 ncssearch = None,
                                                                                 remove_redundance = False,
                                                                                 identity = identity,
                                                                                 ssbridge = False,
                                                                                 peptide_length = 3,
                                                                                 connectivity = False,
                                                                                 weight = 'distance_avg',
                                                                                 strictness_ah = strictness_ah,
                                                                                 strictness_bs = strictness_bs,
                                                                                 list_chains = list_chains,
                                                                                 processed_cvlist=[], interface = True)

                if len(solutions) != 0:
                    print("Interface repeated")
                    to_delete.append(dict_pdbs[perm][dist]['pdbf'])
                    dict_pdbs[perm][dist]['unique'] = 0 #Tag for repeated interfaces
                    # os.remove(dict_pdbs[perm][dist]['pdbf'])
                    continue

            if ind != 0:
                for ind3, perm2 in enumerate(l):
                    # print(perm, perm2, ind3, ind, dist)

                    if ind3 < ind and perm2 in dict_pdbs.keys() and dist in dict_pdbs[perm2].keys():
                        structure, solutions2, chains_dict, lisBigSS = process_structure_fast(
                            pdbsearch=dict_pdbs[perm][dist]['pdbsearch'],
                            pdbf=dict_pdbs[perm2][dist]['pdbf'],
                            pattern=None,
                            cvs_model=dict_pdbs[perm][dist]['cvs_model'],
                            sequence='',
                            ncssearch=None,
                            remove_redundance=False,
                            identity=identity,
                            ssbridge=False,
                            peptide_length=3,
                            connectivity=False,
                            weight='distance_avg',
                            strictness_ah=strictness_ah,
                            strictness_bs=strictness_bs,
                            list_chains=list_chains,
                            processed_cvlist=[], interface = True)

                        if len(solutions2) != 0:
                            print("The interface is repeated with", dict_pdbs[perm][dist]['pdbf'])
                            to_delete.append(dict_pdbs[perm][dist]['pdbf'])
                            dict_pdbs[perm][dist]['unique'] = 0 #Tag for repeated interfaces
                            break
                            # continue
                        else:
                            print("The interface is unique")
                            dict_pdbs[perm][dist]['unique'] = 1 #Tag for unique interfaces
            elif dict_pdbs[perm][dist]['unique'] == 3:
                dict_pdbs[perm][dist]['unique'] = 1  # Tag for unique interfaces

    for file_del in to_delete: os.remove(file_del)
    graph_interfaces = generate_graph_interfaces (dict_pdbs)
    # write_from_graph(graph_interfaces, structure_reference, dir_interfaces, 'interfaces', write_pdb=False,
    #                  write_graphml=True) #TODO: poner bien los paths
    pickle_file = open(os.path.join(dir_interfaces, 'graph_interfaces.pickle'), "wb")
    pickle.dump(graph_interfaces, pickle_file)
    pickle_file.close()

def generate_matrices(model, length_seq, norm_dist=None):

    # Compute CVs and matrices of distances/angles
    strucc, cvs_list = parse_pdb_and_generate_cvs(model)
    F, D, T = decompose_pdb_in_cvs_mat_geom(cvs_list[0], euclidean=True, normalise=True, norm_dist=norm_dist)

    # Dictionary that keep the conversion of CVs id to the residue position
    dict_conversion = {}
    list_res = Bioinformatics.get_list_of_residues(strucc)
    n = 0
    for ires, res in enumerate(list_res):
        try:
            res_num = list_res[ires].get_id()[1]
            res_num2 = list_res[ires + 1].get_id()[1]
            res_num3 = list_res[ires + 2].get_id()[1]
        except:
            break
        if res_num2 == res_num + 1 and res_num3 == res_num + 2:
            dict_conversion[res.get_id()[1]] = n
            n += 1

    # Generate a new matrix for angles/distances. In case the structure misses some residue/CV,
    # the correspondent cell will be filled with Nan value
    new_D = reshape_matrix(dict_conversion, length_seq, D)
    new_T = reshape_matrix(dict_conversion, length_seq, T)

    # #Info: furthest CV in the structure
    # max = numpy.nanmax(new_D)
    # result = numpy.where(new_D == max)
    # print("The furthest CVs are those that start at the residue %i and %i" % (result[0][0] + 1, result[0][1] + 1))

    return new_D, new_T


def reshape_matrix(dict_conversion, length_seq, old_matrix):

    length_seq = length_seq - 2 #The two last CVs cannot be computed

    new_matrix = numpy.empty((length_seq, length_seq))
    new_matrix[:] = numpy.NaN
    for irow, row in enumerate(new_matrix):
        for icolumn, column in enumerate(row):
            if irow + 1 in dict_conversion and icolumn + 1 in dict_conversion:
                try:
                    new_matrix[irow, icolumn] = old_matrix[dict_conversion[irow + 1]][dict_conversion[icolumn + 1]]
                except:
                    new_matrix[irow, icolumn] = numpy.NaN
    return new_matrix

def distance_frobenius(Mr, Mt):

    diff = Mr-Mt

    diff_clean = diff[~numpy.isnan(diff)] #Matrix without Nan values
    core = int(numpy.sqrt(len(diff_clean))) #Core: number of CVs compared

    diff[numpy.isnan(diff)] = 0 #Substitute Nan values for 0 values, which do not influence in the calculations
    diff_transpose = diff.transpose() #Computing transpose matrix
    product = numpy.dot(diff, diff_transpose)
    trace = product.trace()
    dist = numpy.sqrt(trace)

    return dist, core


def saving_results_df(dict_df, reference, dist_norm, dist_frob):

    reference = os.path.basename(reference)[:-4]

    if not reference in dict_df[dist_norm]:
        dict_df[dist_norm][reference] = []
    dict_df[dist_norm][reference].append(dist_frob)

    return dict_df


def generate_heatmap_plot(matrix, title, name, vmin=0, vmax=None):
    color_map = plt.cm.get_cmap('RdYlGn')
    reversed_color_map = color_map.reversed()

    fig = plt.figure()
    fig.suptitle(title, fontsize=20)
    plt.imshow(matrix, cmap=reversed_color_map, interpolation='nearest', vmin=vmin, vmax=vmax)
    plt.colorbar()

    plt.savefig((name + '.png'), bbox_inches='tight', dpi=150)
    return (os.path.join(os.getcwd(), "plots", name + '.png'))


def generate_scatter_plot(list_x, list_y, name):

    fig = plt.figure()
    fig.suptitle(os.path.basename(name), fontsize=20)

    plt.plot(list_x, list_y)
    plt.xticks(list_x, list_x)
    plt.xlabel('Distance normalisation (Å)')
    plt.ylabel('Frobenius distance')

    plt.savefig(os.path.join(name + '.png'), dpi=150)


def plot_dendrogram(model, targets, **kwargs):

    fig = plt.figure()
    plt.title("Hierarchical Clustering Dendrogram")

    # Create linkage matrix and then plot the dendrogram
    # Create the counts of samples under each node
    counts = numpy.zeros(model.children_.shape[0])
    n_samples = len(model.labels_)

    for i, merge in enumerate(model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = numpy.column_stack([model.children_, model.distances_, counts]).astype(float)

    # Plot the corresponding dendrogram
    scipy.cluster.hierarchy.dendrogram(linkage_matrix, **kwargs, labels=[target for target in targets])
    plt.xlabel("Number of points in node (or index of point if no parenthesis).")
    plt.savefig(os.path.join('clustering_dendogram_angles.png'), dpi=150)

    fl = scipy.cluster.hierarchy.fcluster(linkage_matrix, max(model.distances_)-0.1, criterion='distance', depth=190)
    return_dict = {}
    for i, element in enumerate(fl):
        return_dict.setdefault(element, []).append(list(targets.keys())[i])
    return list(return_dict.values())

def frobenius(references, targets, norm_dist=None, min_dist=None, max_dist=None, step_dist=None,
              write_plot=False, write_matrix=False):

    wd=os.getcwd()
    references.sort()
    targets.sort()

    #These list will save the results for returning when function is used as a library
    list_references = []
    list_targets =[]
    list_distances_norm = []
    list_frobenius_angles = []
    list_frobenius_distances = []
    list_plot_ang = []
    list_plot_dist = []
    list_core = []

    #Directory for saving matrices
    if write_matrix:
        if os.path.exists('./matrices'): shutil.rmtree('./matrices')
        os.makedirs('./matrices')

    # Distance normalisation
    if norm_dist:  # In case the user provides a value for distance normalisation
        distances = [norm_dist]
    elif min_dist and max_dist and step_dist:  # In case the user provide a range
        distances = [x for x in range(min_dist, max_dist + 1, step_dist)]
    else: #Default will be 30Å if distance for normalisation is not provided
        print("\nINFO: A distance of 30A will be taken for normalisation\n")
        distances = [30]

    #Dictionaries for saving results for plots/dataframes/dendogram
    dict_norm = {}
    dict_df_ang = {dist:{} for dist in distances}
    dict_df_dist = {dist:{} for dist in distances}
    dict_matrix_ang = {}


    for reference in references:
        #Getting the last residue of the structure, ignoring HETATM, HOH, etc
        #It will be used for setting the matrix dimension
        strucc = Bioinformatics.get_structure("stru", reference)
        chains = Bio.PDB.Selection.unfold_entities(strucc, "C")
        for residue in Bio.PDB.Selection.unfold_entities(chains[0], "R"):
            if (residue.get_resname().upper() in Bioinformatics.AAList) and (residue.has_id("CA")) and (
                    residue.has_id("O")) and (residue.has_id("C")) and (residue.has_id("N")):
                length_seq = residue.get_id()[1]

        #Directory for saving plots
        if write_plot:
            dir_out = os.path.join(wd, os.path.basename(reference)[:-4] + '_ref')
            if os.path.exists(dir_out): shutil.rmtree(dir_out)
            os.mkdir(dir_out)

        for i, dist in enumerate(distances):
            print("---------------------------------------------------------------------")
            print("Normalising distances at %s A: \n" % (dist))
            print("Analysing reference file: ", os.path.basename(reference))
            Dr, Tr = generate_matrices(model=reference, length_seq=length_seq, norm_dist=dist)

            # if write_plot:
                # generate_heatmap_plot(matrix=Tr, name= os.path.join(dir_out, 'ang_' + os.path.basename(reference)[:-4]))
                # generate_heatmap_plot(matrix=Dr, name= os.path.join(dir_out, 'dist_' + os.path.basename(reference)[:-4]))

            if write_matrix:
                numpy.save(os.path.join(wd, 'matrices', os.path.basename(reference)[:-4] + '_dist_' + str(dist) + 'A'), Dr)
                numpy.save(os.path.join(wd, 'matrices', os.path.basename(reference)[:-4] + '_ang'), Tr)

            for target in targets:
                if target.endswith('.pdb'):
                    print("\nAnalysing target file: ", os.path.basename(target))
                    Dt, Tt = generate_matrices(model=target, length_seq=length_seq, norm_dist=dist)

                    # if write_plot:
                    #     generate_heatmap_plot(matrix=Tt, name=os.path.join(dir_out, 'ang_' + os.path.basename(target)[:-4]))
                    #     generate_heatmap_plot(matrix=Dt, name=os.path.join(dir_out,'dist_' + os.path.basename(target)[:-4]))

                    if write_matrix:
                        numpy.save(os.path.join(wd, 'matrices',os.path.basename(target)[:-4] + '_dist_' + str(dist) + 'A'), Dt)
                        numpy.save(os.path.join(wd, 'matrices',os.path.basename(target)[:-4] + '_ang'), Tt)

                    #Computing differences between matrices form reference and target
                    D_difference = numpy.absolute(Dr - Dt)
                    T_difference = numpy.absolute(Tr - Tt)

                    if write_plot:
                        name_diff = os.path.join(dir_out, os.path.basename(reference[:-4] + '_VS_' + os.path.basename(target)[:-4]))
                        title = f'DIST: {os.path.basename(reference)[:-4]} VS {os.path.basename(target)[:-4]}'
                        path_plot_dist = generate_heatmap_plot(matrix=D_difference, title=title, name=name_diff + '_dist_' + str(dist) + 'A' , vmin=0, vmax=1)
                        list_plot_dist.append(path_plot_dist)
                        title = f'ANG: {os.path.basename(reference)[:-4]} VS {os.path.basename(target)[:-4]}'
                        path_plot_ang = generate_heatmap_plot(matrix=T_difference, title=title, name=name_diff + '_ang' , vmin=0, vmax=1)
                        list_plot_ang.append(path_plot_ang)

                    #Computing Frobenius distance
                    dist_frobenius_ang, core = distance_frobenius(Tr, Tt)
                    dist_frobenius_dist, core = distance_frobenius(Dr, Dt)
                    print("    Distance Frobenius for CV angles is: ",dist_frobenius_ang.round(2))
                    print("    Distance Frobenius for CV distances is: ", dist_frobenius_dist.round(2))
                    print("    Core is : %s CV"%(core))

                    #Saving the results in dictionary for plotting as a dataframe
                    dict_df_ang = saving_results_df(dict_df_ang, reference, dist, dist_frobenius_ang.round(2))
                    dict_df_dist = saving_results_df(dict_df_dist, reference, dist, dist_frobenius_dist.round(2))

                    #Saving results in lists for returning when function is used as a library
                    list_references.append(os.path.basename(reference)[:-4])
                    list_targets.append(os.path.basename(target)[:-4])
                    list_distances_norm.append(dist)
                    list_frobenius_angles.append(dist_frobenius_ang.round(2))
                    list_frobenius_distances.append(dist_frobenius_dist.round(2))
                    list_core.append(core)

                    if not os.path.basename(target)[:-4] in dict_norm:
                        dict_norm[os.path.basename(target)[:-4]] = [0 for x in distances]
                    dict_norm[os.path.basename(target)[:-4]][i] = dist_frobenius_dist.round(2)

                    if not os.path.basename(target)[:-4] in dict_matrix_ang:
                        new_Tt = Tt.flatten()
                        new_Tt[numpy.isnan(new_Tt)] = 0
                        dict_matrix_ang[os.path.basename(target)[:-4]]= new_Tt

            print('\n')

        #Writting plot for distance normalisation
        if len(distances) > 1 and write_plot:
            for model in dict_norm:
                generate_scatter_plot(list_x=distances, list_y=dict_norm[model], name=os.path.join(dir_out, os.path.basename(reference)[:-4] + '_VS_' + model + '_dist_norm'))


    # Printing results as a dataframe and saving them as .csv file
    for dist in dict_df_dist:
        df_angles = pandas.DataFrame({})
        df_distances = pandas.DataFrame({})

        for ref in dict_df_dist[dist]:
            df_angles[ref] = dict_df_ang[dist][ref]
            df_distances[ref] = dict_df_dist[dist][ref]

        df_angles.index = [os.path.basename(x)[:-4] for x in targets]
        df_distances.index = [os.path.basename(x)[:-4] for x in targets]

        print('\n\nFrobenius distances of CV distances (%s A):' % (dist))
        print(df_distances)
        df_distances.to_csv('frobenius_distances_' + str(dist) + '.csv', sep='\t')

    print('\n\nFrobenius distances of CV angles:')
    print(df_angles)
    df_angles.to_csv('frobenius_angles.csv', sep='\t')
    print('\n\n')

    # Clustering and writtting dendogram
    list_matrix = [dict_matrix_ang[matrix] for matrix in dict_matrix_ang]

    dendogram_list = []
    if len(list_matrix) > 1:
        model = AgglomerativeClustering(distance_threshold=0, n_clusters=None) #setting distance_threshold=0 ensures we compute the full tree
        model = model.fit(list_matrix)
        dendogram_list = plot_dendrogram(model, truncate_mode="level", p=2, targets=dict_matrix_ang) #plot the top three levels of the dendrogram
    else:
        print('\n\nNot possible to calculate the dendrogram with just one sample')

    return list_references, list_targets, list_core, list_distances_norm, list_frobenius_angles, list_frobenius_distances, list_plot_ang, list_plot_dist, dendogram_list


def frobenius_submatrices(path_ref, path_tar, residues_ref=None, residues_tar=None, write_plot=False, title=''):

    print('Reading matrix files: ', os.path.basename(path_ref), os.path.basename(path_tar))
    matrix_ref = numpy.load(path_ref)
    matrix_tar = numpy.load(path_tar)


    if residues_ref:
        matrix_ref = select_submatrix(matrix_ref, residues_ref)
        print('Reference submatrix will be computed for residues %s' % (residues_ref))

    if residues_tar:
        matrix_tar = select_submatrix(matrix_tar, residues_tar)
        print('Target submatrix will be computed for residues:\n    %s' % (residues_tar))

    difference = numpy.absolute(matrix_ref - matrix_tar)

    path_plot = ''
    if write_plot:
        if not os.path.exists('./plots'): os.makedirs('./plots')
        name_diff = os.path.basename(path_ref)[:-4] + '_VS_' + os.path.basename(path_tar)[:-4] + '_cut'
        if title == '':
            title = f'CUT: {os.path.basename(path_ref)[:-4]} VS {os.path.basename(path_tar)[:-4]}'
        path_plot = generate_heatmap_plot(matrix=difference, title=title, name=name_diff, vmin=0, vmax=1)


    dist_frobenius, core = distance_frobenius(matrix_ref, matrix_tar)
    print("\n    Distance Frobenius is: ", dist_frobenius.round(2))
    print("    Core is : %s CV\n" % (core))

    return ([dist_frobenius.round(2), core, path_plot])

def select_submatrix(matrix, residues_tar):

    new_matrix = numpy.empty(matrix.shape)
    new_matrix[:] = numpy.NaN

    for irow, row in enumerate(new_matrix):
        for icolumn, column in enumerate(row):
            if (set([irow + 1, irow+2, irow+3]).issubset(set(residues_tar))) and (set([icolumn + 1, icolumn+2, icolumn+3]).issubset(set(residues_tar))): #Correspondance between CVs and residues
                try:
                    new_matrix[irow, icolumn] = matrix[irow][icolumn]
                except:
                    new_matrix[irow, icolumn] = numpy.NaN

    return(new_matrix)


def remove_unconnected_fragments(reference, strictness_ah, strictness_bs, peptide_length, write_pdb=False, domains=False,
                                  write_graphml=False, write_graph_img=False, path_base="./"):

    if strictness_ah < 0:
        strictness_ah = 0.45
    if strictness_bs < 0:
        strictness_bs = 0.20

    weight = "distance_avg"
    pep_len = int(peptide_length)

    pdb_search_in = ""
    f = open(reference, "r")
    all_lines = f.readlines()
    f.close()
    for line in all_lines:
        if line[:6] in ["TITLE ", "CRYST1", "SCALE1", "SCALE2", "SCALE3"]:
            pdb_search_in += line

    graph_full_reference, structure_reference, matrix_reference, cvs_reference, highd_reference = annotate_pdb_model_with_aleph(
        reference, weight=weight, strictness_ah=strictness_ah, strictness_bs=strictness_bs, peptide_length=pep_len,
        write_pdb=False, domains=False)
    graph_no_coil_reference = graph_full_reference.vs.select(sstype_in=["ah", "bs"]).subgraph()

    edges_filter = graph_no_coil_reference.es.select(min_le=20)
    graph_no_coil_reference_filtered = graph_no_coil_reference.subgraph_edges(edges_filter, delete_vertices=True)


    list_compact = []
    for frag in graph_no_coil_reference_filtered.vs():
        list_compact.append(frag['reslist'][0])

    list_frag_sorted = []
    for frag in graph_full_reference.vs():
        list_frag_sorted.append([frag['sstype'], list(frag['reslist'][0]), frag.index])

    list_frag_sorted = sorted(list_frag_sorted, key=lambda x: x[1]) #Sort by sequence

    prev = False
    list_vertex_delete = []
    print('\n\n\n')
    for frag in list_frag_sorted:
        if frag[0] == 'coil' and prev == False:
            # print('no', 'coil', frag[1][3][1], prev)
            list_vertex_delete.append(frag[2])

        elif frag[0] == 'coil' and prev == True:
            # print('yes', 'coil', frag[1][3][1], prev)
            continue

        elif frag[0] in ['ah', 'bs'] and frag[1] in list_compact:
            # print('yes', frag[1][3][1], frag[0], prev)
            prev = True
        else:
            # print('no',frag[1][3][1], frag[0], prev)
            prev = False
            list_vertex_delete.append(frag[2])

    graph_full_reference.delete_vertices(list_vertex_delete)
    write_from_graph(graph_full_reference, structure_reference, './', os.path.basename(reference)[:-4] + '_input_search_compact', write_pdb=True, write_graphml=False)


#######################################################################################################
#                                               MAIN                                                  #
#######################################################################################################

def main():
    start_time = time.time()

    head1 = """
                                      .------------------------------------------.
                                      |            _      ______ _____  _    _   |
                                      |      /\   | |    |  ____|  __ \| |  | |  |
                                      |     /  \  | |    | |__  | |__) | |__| |  |
                                      |    / /\ \ | |    |  __| |  ___/|  __  |  |
                                      |   / ____ \| |____| |____| |    | |  | |  |
                                      |  /_/    \_\______|______|_|    |_|  |_|  | 
                                      |                                          | 
                                      #------------------------------------------#
                                                | v. 1.3.1  -- 12/2020 |
                    """
    cprint(head1, 'cyan')
    print("""
        In case this result is helpful, please, cite:

        ALEPH: a network-oriented approach for the generation of fragment-based libraries and for structure interpretation
        Medina A, Triviño J, Borges RJB, Millan C , Usón I and Sammito MD
        (2019) Acta Cryst. D Study Weekend
        
        Exploiting tertiary structure through local folds for ab initio phasing
        Sammito M, Millán C, Rodríguez DD, M. de Ilarduya I, Meindl K, De Marino I, Petrillo G, Buey RM, de Pereda JM,
        Zeth K, Sheldrick GM and Usón I
        (2013) Nat Methods. 10, 1099-1101.
        """)
    print("Email support: ", colored("alephbugs@gmail.com", 'blue'))

    # List of arguments accepted in the command line
    parser = argparse.ArgumentParser(description='Command line options for ALEPH')
    general_group = argparse.ArgumentParser(add_help=False)
    subparsers = parser.add_subparsers(title='List of features in ALEPH', description='Functions implemented in ALEPH', help='Help for each feature')

    general_group.add_argument("--width_pic", help="Width in inches for pictures. Default 100.0", type=float, default=100.0)
    general_group.add_argument("--height_pic", help="Height in inches for pictures. Default 20.0", type=float, default=20.0)
    general_group.add_argument("--pack_beta_sheet",
                        help="Do not break a beta sheet in two community clusters unless is really necessary. Default: False",
                        action='store_true', default=False)
    general_group.add_argument("--homogeneity",
                        help="Favourite homogeneous clusters. Default: False",
                        action='store_true', default=False)
    general_group.add_argument("--max_ah_dist", help="Maximum distance allowed among ah cvs in the graph. Default 20.0",
                        type=float, default=20.0)
    general_group.add_argument("--min_ah_dist", help="Minimum distance allowed among ah cvs in the graph. Default  0.0",
                        type=float, default=0.0)
    general_group.add_argument("--max_bs_dist", help="Maximum distance allowed among bs cvs in the graph. Default 15.0",
                        type=float, default=15.0)
    general_group.add_argument("--min_bs_dist", help="Minimum distance allowed among bs cvs in the graph. Default  0.0",
                        type=float, default=0.0)
    general_group.add_argument("--rmsd_thresh", help="Rmsd threshold to accept a superposition. Default: 1.5",
                        type=float, default=1.5)
    general_group.add_argument("--write_graphml", action="store_true", help="Write graphml files. Default: False",
                        default=False)
    general_group.add_argument("--strictness_ah", help="strictness parameter threshold for accepting ah CVs. Default: 0.50",
                        type=float, default=0.50)
    general_group.add_argument("--strictness_bs", help="strictness parameter threshold for accepting bs CVs. Default: 0.30",
                        type=float, default=0.30)
    general_group.add_argument("--peptide_length", help="Define the peptide length for computing a CVs. Default: 3", type=int,
                        default=3)

    parser_cluster = subparsers.add_parser("cluster", help='Cluster a library from ALEPH', parents=[general_group])
    parser_cluster.add_argument("--directory_database", help="Directory with pdb file models representing the library", action='store', required=True)
    parser_cluster.add_argument("--rmsd_clustering",
                                        help="Rmsd threshold for geometrical clustering of the library. Default: 1.5",
                                        type=float, default=1.5)
    parser_cluster.add_argument("--clustering_mode", help="Clustering algorithm. Default: no_clustering",
                                        type=str, default="no_clustering",
                                        choices=["rmsd", "rmsd_range", "random_sampling", "no_clustering"])
    parser_cluster.add_argument("--number_of_ranges",
                                        help="If rmsd_range clustering algorithm is activated, it specifies the number of ranges to group the extracted models. Default: 500",
                                        type=int, default=500)
    parser_cluster.add_argument("--number_of_clusters",
                                        help="If the rmsd_range or random sampling algorithm is activated, it specifies the absolute number of representative models extracted from the library. Default: 7000",
                                        type=int, default=7000)
    parser_cluster.add_argument("--exclude_residues_superpose",
                                        help="Number of residues to possibly exclude from the superposition core. Default: 0 (== No exclusion)",
                                        type=int, default=0)
    parser_cluster.add_argument("--ssbridge", action="store_true",
                                        help="Check for disulphide bridges. Default: False", default=False)
    parser_cluster.add_argument("--nilges", help="Cycles of iteration for the nilges algorithm. Default: 10",
                                        type=int, default=10)
    parser_cluster.set_defaults(func=perform_cluster_starter)

    parser_superpose = subparsers.add_parser("superpose", help='Superpose two pdb models', parents=[general_group])
    parser_superpose.add_argument("--reference", help="Input reference pdb model (that is the complete structure)", required=True)
    parser_superpose.add_argument("--target", help="Input target pdb model (that is the small fold)", required=False)
    parser_superpose.add_argument("--targets", help="Input target pdb model (that are small folds)", required=False)
    parser_superpose.add_argument("--nilges", help="Cycles of iteration for the nilges algorithm. Default: 10", type=int, default=10)
    parser_superpose.add_argument("-C", "--score_intra_fragment", help="Global geometrical match between template and extracted individual fragments expressed as score percentage.", type=int, required=False)
    parser_superpose.add_argument("-J", "--score_inter_fragments", help="Global geometrical match between template and extracted fold expressed as score percentage.", type=int, required=False)
    parser_superpose.add_argument("--reverse", action="store_true", help="Superpose the reference on top of the target(s). Default: False", default=False)
    parser_superpose.add_argument("--criterium_selection_core", type=str, help="What is the criteria to be used to compute the core_percentage. Default: residues", default="residues", choices=["residues","secondary_structures"])
    parser_superpose.add_argument("--multimer", action="store_false", help="Remove chain redundancy. Default: True", default=True)
    parser_superpose.add_argument("--similarity_intra_chain", help="Global geometrical match between template fragments and extracted fragments expressed as a number between 1 (hardly similar) "
                                                                         "and 5 (almost identical). Default: 3", type=int, default=3)
    parser_superpose.add_argument("--similarity_inter_chains", help="Global geometrical match between interactions observed among template fragments and extracted fold expressed as a number between 1 (hardly similar)"
                                                                         " and 5 (almost identical). Default: 3", type=int, default=3)
    parser_superpose.add_argument("--representative", action="store_true", help="Superpose each target to reference only one time, the one with the lowest rmsd. Default: False", default=False)
    parser_superpose.add_argument("--interface", action="store_true", help="Superpose each target also to reference symmetry equivalents. Default: False", default=False)
    parser_superpose.add_argument("--sidechains", action="store_true", help="Output models with side chains. Default: False", default=False)
    parser_superpose.add_argument("--verbose", action="store_true", help="Verbose output. Default: False",
                                        default=False)
    parser_superpose.set_defaults(func=perform_superposition_starter)

    parser_annotate = subparsers.add_parser("annotate", help="Annotate pdbmodel with CVs and produces text and image reports", parents=[general_group])
    parser_annotate.add_argument("--pdbmodel", help="Input a pdb model ", required=True)
    parser_annotate.add_argument("--hhr_file", help="Input the .hhr file containing the multiple alignment", type=str, default='')
    parser_annotate.add_argument("--verbose", action="store_true", help="Verbose output. Default: False", default=False)
    parser_annotate.add_argument("--compactness", action="store_true", help="Annotate compact fragments. Default: False", default=False)
    parser_annotate.set_defaults(func=annotate_pdb_model_starter)

    parser_decompose = subparsers.add_parser("decompose", help="Compute community clustering for decomposition in structural units", parents=[general_group])
    parser_decompose.add_argument("--pdbmodel", help="Input a pdb model ", required=True)
    parser_decompose.add_argument("--algorithm", help="Algorithm for the community clustering procedure.", type=str, default="fastgreedy", choices=['fastgreedy', 'infomap', 'eigenvectors', 'label_propagation','community_multilevel', 'edge_betweenness', 'spinglass', 'walktrap'])
    parser_decompose.add_argument("--hhr_file", help="Input the .hhr file containing the multiple alignment", type=str, default='')
    parser_decompose.add_argument("--write_graph_img", action="store_true", help="Write image for graphs. Default: False", default=False)
    parser_decompose.set_defaults(func=decompose_by_community_clustering_starter)

    parser_findfolds = subparsers.add_parser("find_folds", help="Search and exctracts folds in a protein structure", parents=[general_group])
    parser_findfolds.add_argument("--pdbmodel", help="Input a pdb model ", required=True)
    parser_findfolds.add_argument("--domains", action="store_true", help="Decomposition intending to reveal domains. Default: False", default=False)
    parser_findfolds.add_argument("--compactness", action="store_true", help="Annotate compact fragments. Default: False", default=False)
    parser_findfolds.add_argument("--write_graph_img", action="store_true", help="Write image for graphs. Default: False", default=False)
    parser_findfolds.set_defaults(func=find_local_folds_in_the_graph_starter)

    parser_generatelibrary= subparsers.add_parser("generate_library", help="Generate a library of the given fold superposed to the template ready for being used in ARCIMBOLDO_BORGES", parents=[general_group])
    parser_generatelibrary.add_argument("--pdbmodel", help="Input a pdb model ")
    parser_generatelibrary.add_argument("--use_model_as_it_is", action="store_true", help="Verbose output. Default: False", default=False)

    parser_generatelibrary.add_argument("--targz",    help="Read all input from a tar.gz pre-formatted with ALEPH for a Grid.")

    parser_generatelibrary.add_argument("--directory_database", help="Directory with pdb file of the deposited structures", action='store', required=False)
    parser_generatelibrary.add_argument("--cath_id", help="Extract a database from the given cath id", required=False)
    parser_generatelibrary.add_argument("--target_sequence", help="Extract a database from the given target sequence", required=False)
    parser_generatelibrary.add_argument("--exclude_sequence", help="Avoid the extraction of models from any chain that align with this sequence with a s.i. >= 90%%", required=False)
    parallelization = parser_generatelibrary.add_mutually_exclusive_group()
    parallelization.add_argument("--supercomputer", help="Nodefile for the supercomputer")
    parallelization.add_argument("--local_grid", help="Path to the setup.bor to start jobs in local grid")
    parallelization.add_argument("--remote_grid", help="Path to the setup.bor to start jobs in remote grid")

    parser_generatelibrary.add_argument("--work_directory", help="Working Directory path. Default is ./", default="./")
    parser_generatelibrary.add_argument("--criterium_selection_core", type=str, help="What is the criteria to be used to compute the core_percentage. Default: residues", default="residues", choices=["residues","secondary_structures"])
    parser_generatelibrary.add_argument("-C", "--score_intra_fragment", help="Global geometrical match between template and extracted individual fragments expressed as score percentage.", type=int, required=False)
    parser_generatelibrary.add_argument("--c_angle", help="Percentage of agreement for internal angles in a fragment. Default: 95", type=int, default=-1)
    parser_generatelibrary.add_argument("--c_dist", help="Percentage of agreement for internal distances in a fragment. Default: 95", type=int, default=-1)
    parser_generatelibrary.add_argument("--c_angle_dist", help="Percentage of agreement for internal distance angles in a fragment. Default: 95", type=int, default=-1)
    parser_generatelibrary.add_argument("--c_cvl_diff", help="Percentage of agreement for internal CVL differences observed in a fragment. Default: 95", type=int, default=-1)
    parser_generatelibrary.add_argument("-J", "--score_inter_fragments", help="Global geometrical match between template and extracted fold expressed as score percentage",type=int, required=False)
    parser_generatelibrary.add_argument("--j_angle", help="Percentage of agreement for external angles between two fragments. Default: 90", type=int, default=-1)
    parser_generatelibrary.add_argument("--j_dist", help="Percentage of agreement for external distances between two fragments. Default: 90", type=int, default=-1)
    parser_generatelibrary.add_argument("--j_angle_dist", help="Percentage of agreement for external distance angles between two fragments. Default: 90", type=int, default=-1)
    parser_generatelibrary.add_argument("--j_cvl_diff", help="Percentage of agreement for external CVL differences observed between two fragments. Default: 90", type=int, default=-1)
    parser_generatelibrary.add_argument("--similarity_intra_chain", help="Global geometrical match between template fragments and extracted fragments expressed as a number between 1 (hardly similar) "
                                                                         "and 5 (almost identical). Default: 3", type=int, default=3)
    parser_generatelibrary.add_argument("--similarity_inter_chains", help="Global geometrical match between interactions observed among template fragments and extracted fold expressed as a number between 1 (hardly similar)"
                                                                         " and 5 (almost identical). Default: 3", type=int, default=3)


    parser_generatelibrary.add_argument("--verbose", action="store_true", help="Verbose output. Default: False", default=False)
    parser_generatelibrary.add_argument("--sidechains", action="store_true", help="Output models with side chains. Default: False", default=False)
    parser_generatelibrary.add_argument("--sequence", help="Require a specific sequence in the output model to match the template. Complete template sequence needs to be given; X marks unspecified residues.", type=str, default="")
    parser_generatelibrary.add_argument("--ncssearch", action="store_true", help="Extract local folds also from NCS relative copies. Default: False", default=False)
    parser_generatelibrary.add_argument("--multimer", action="store_false", help="Remove chain redundancy unless NCS is set. Default: True", default=True)
    parser_generatelibrary.add_argument("--force_core", help="Number of parallel processes. Default: -1 (== #cores machine)", type=int, default=-1)
    parser_generatelibrary.add_argument("--rmsd_min",  help="Minimum rmsd against the template. Default: 0.0 (== extract identical)", type=float, default=0.0)
    parser_generatelibrary.add_argument("--rmsd_max",  help="Maximum rmsd against the template. Default: 6.0 (== model difference no grater then 6.0 A)", type=float, default=6.0)
    parser_generatelibrary.add_argument("--rmsd_clustering", help="Rmsd threshold for geometrical clustering of the library. Default: 1.5", type=float, default=1.5)
    parser_generatelibrary.add_argument("--clustering_mode", help="Clustering algorithm. Default: no_clustering", type=str, default="no_clustering", choices=["rmsd", "rmsd_range", "random_sampling", "no_clustering"])
    parser_generatelibrary.add_argument("--number_of_ranges", help="If rmsd_range clustering algorithm is activated, it specifies the number of ranges to group the extracted models. Default: 500", type=int, default=500)
    parser_generatelibrary.add_argument("--number_of_clusters",help="If the rmsd_range or random sampling algorithm is activated, it specifies the absolute number of representative models extracted from the library. Default: 7000", type=int, default=7000)
    parser_generatelibrary.add_argument("--exclude_residues_superpose", help="Number of residues to possibly exclude from the superposition core. Default: 0 (== No exclusion)", type=int, default=0)
    parser_generatelibrary.add_argument("--ssbridge", action="store_true", help="Check for disulphide bridges. Default: False", default=False)
    parser_generatelibrary.add_argument("--nilges", help="Cycles of iteration for the nilges algorithm. Default: 10", type=int, default=10)
    parser_generatelibrary.add_argument("--enhance_fold", action="store_true",help="Use the minimum fragment length in the template as fragment size for computing CVLs",  default=False)
    parser_generatelibrary.add_argument("--representative", action="store_true", help="For each structure in the PDB database extracts only one model, the one with the lowest rmsd. Default: False",  default=False)
    parser_generatelibrary.add_argument("--remove_coil", action="store_true",  help="Remove coil regions from template before searching. Default: False", default=False)
    parser_generatelibrary.add_argument("--connectivity", action="store_true", help="Fragments extracted must have the same sequence order as template.  Default: False", default=False)
    parser_generatelibrary.add_argument("--test", action="store_true", help="Test with a reduced sample of models to check parameterisation.  Default: False", default=False)
    parser_generatelibrary.add_argument("--interface", action="store_true", help="Extract local folds also from the symmetry equivalents. Default: False", default=False)

    parser_generatelibrary.set_defaults(func=generate_library_starter)

    parser_interfaces = subparsers.add_parser("analyse_interfaces", help="Generate interfaces from template in a given radius with spack and analyse them",parents=[general_group])
    parser_interfaces.add_argument("--pdbmodel", help="Input a pdb model ", required=True)
    parser_interfaces.add_argument("--min_distance", help="Minimum radius in angstrom for generating the symmetry copies. Default: None", type=int, default=None)
    parser_interfaces.add_argument("--max_distance", help="Maximum radius in angstrom for generating the symmetry copies. Default: None", type=int, default=None)
    parser_interfaces.add_argument("--edges_distance", help="Distance in angstrom for filtering the edges of the graph. Default: 20", type=int, default=20)
    parser_interfaces.add_argument("--verbose", action="store_true", help="Verbose output. Default: False", default=False)
    parser_interfaces.set_defaults(func=analyse_interfaces_starter)

    parser_frobenius = subparsers.add_parser("frobenius", help="Analyse geometry and compute Frobenius distance between models")
    parser_frobenius.add_argument("--reference", help="Input reference pdb model", required=False)
    parser_frobenius.add_argument("--references", help="Path to the directory containing the reference pdb model", required=False)
    parser_frobenius.add_argument("--targets", help="Path to the directory containing target pdb file(s)", required=True)
    parser_frobenius.add_argument("--distance", help="Value for normalisation of the distance", type=int, required=False)
    parser_frobenius.add_argument("--distance_max", help="Maximum value for normalisation of the distance", type=int, required=False)
    parser_frobenius.add_argument("--distance_min", help="Maximum value for normalisation of the distance", type=int, required=False)
    parser_frobenius.add_argument("--distance_step", help="Step value for normalisation of the distance", type=int, required=False)
    parser_frobenius.add_argument("--write_plot", action="store_true",  help="Write plots for difference between matrices", default=False, required=False)
    parser_frobenius.add_argument("--write_matrix", action="store_true", help="Write .pkl files containing matrices", default=False, required=False)
    parser_frobenius.set_defaults(func=frobenius_starter)

    parser_compactness = subparsers.add_parser("annotate_compactness", help="Remove unnconnected regions", parents=[general_group])
    parser_compactness.add_argument("--pdbmodel", help="Input a pdb model ")
    parser_compactness.set_defaults(func=remove_unconnected_fragments_starter)

    print("*******************************************COMMAND LINE**************************************************")
    print(" ".join(sys.argv))
    print("*********************************************************************************************************")

    args = parser.parse_args()
    try:
        args.func(args)
    except AttributeError:
        parser.print_usage()

    print("Time elapsed: {:.2f}s".format(time.time() - start_time))

    t = datetime.datetime.now()
    epcSec = time.mktime(t.timetuple())
    now = datetime.datetime.fromtimestamp(epcSec)
    print("" + str(now.ctime()))

    print("Job ended with success")

if __name__ == "__main__":
    try:
        main()
    except:
        sys.exit(1)
