import os, json, sys


sys.path.append('..')







from src.bioiain.aleph import CVMatrix, FragmentedStructure
from src.bioiain.base import BIStructure



from src.bioiain.biopython.imports import *
pdb = "6F63"
#pdb = "5LXN"
#pdb = "1M2Z"

pdb_folder = downloadPDB("./data", pdb_list=[pdb], list_name="ccs", file_format="cif")


log("start", "aleph.py")

from src.bioiain.aleph import *
from src.bioiain.biopython import imports
log("header", pdb)

path = os.path.join(pdb_folder, f"{pdb}.cif")
structure = FragmentedStructure.from_file(path)
log(1, structure)
structure.fragment_with_aleph()
cvmap = structure.map_cvectors()
cvmap.save_fig()
cvmap.calculate_neighbours()


exit()
graph, strucc, matrix, cvs_list, highd = core.ALEPH.annotate_pdb_model_with_aleph(pdb_model=path,
            weight="distance_avg",
            strictness_ah=0.45,
            strictness_bs=0.20,
            peptide_length=3,
            write_pdb=False)



print("graph", type(graph), graph.vs)
for v in graph.vs:
    print(v["cvids"])
    #print(dir(v))
    #exit()
print("strucc", type(strucc), len(strucc))
print("matrix", type(matrix), len(matrix))
#print(matrix.keys())

print("cvs_list", type(cvs_list), len(cvs_list))
print("highd", type(highd), highd)


exit()
pdb_folder = downloadPDB("./data", pdb_list=[pdb], list_name="ccs")

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
[print(c.closest)for c in cvmatrix.vectors]












