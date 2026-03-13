import os, json, sys
sys.path.append('..')

from src.bioiain.aleph import CVMatrix

from src.bioiain.biopython.imports import *

pdb_folder = downloadPDB("./data", pdb_list=["5LXN"], list_name="ccs")
structure = loadPDB(os.path.join(pdb_folder, "5LXN.cif"))

print(structure)


for chain in structure.get_chains():
    print(chain)

    #print(chain.residues())

    cvectors = chain._calculate_cvectors()

    cvmatrix = CVMatrix(cvectors) 

    cvmatrix.show("d")
    cvmatrix.show("a")
    cvmatrix.show("t1")
    cvmatrix.show("t2")








