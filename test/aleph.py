import os, json, sys
sys.path.append('..')

from src.bioiain.aleph import *

from src.bioiain.biopython.imports import *

pdb_folder = downloadPDB("./data", pdb_list=["5LXN"], list_name="ccs")
structure = loadPDB(os.path.join(pdb_folder, "5LXN.cif"))

print(structure)


for chain in structure.get_chains():
    print(chain)

    #print(chain.residues())

    chain._calculate_cvectors()








