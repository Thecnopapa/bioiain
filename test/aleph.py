import os, json, sys
sys.path.append('..')

from src.bioiain.aleph import CVMatrix

from src.bioiain.biopython.imports import *
pdb = "6F63"
#pdb = "5LXN"
pdb_folder = downloadPDB("./data", pdb_list=[pdb], list_name="ccs")
structure = loadPDB(os.path.join(pdb_folder, f"{pdb}.cif"))

print(structure)


cvectors = []

for chain in structure.get_chains():
    print(chain)

    #print(chain.residues())

    cvectors.extend(chain.cvectors())

cvmatrix = CVMatrix(cvectors) 


fig_folder = f"./cvmatrixes/{pdb}"

cvmatrix.save_fig("d", save_folder=fig_folder)
cvmatrix.save_fig("a", save_folder=fig_folder)
cvmatrix.save_fig("t1", save_folder=fig_folder)
cvmatrix.save_fig("t2", save_folder=fig_folder)


cvmatrix.calculate_neighbours()
print(c.closest for c in cvmatrix.vectors)








