

import os
import sys

sys.path.append('..')

from src.bioiain.biopython.imports import *
from src.bioiain.symmetries.crystal import *
from src.bioiain.tools.DSSP import DSSP


file_folder = downloadPDB("./data",
                          list_name="mari-weird",
                          file_format="cif",
                          pdb_list=["1M2Z"],
                          overwrite=False)

print(file_folder)
print(os.listdir(file_folder))


dssp = DSSP(dssp_cmd="mkdssp")
print(dssp)
for  file in os.listdir(file_folder):
    if "1M2Z" not in file:
        continue
    monomers = get_monomers(file, file_folder)
    print(monomers)
    for monomer in monomers:
        #[print(a) for a in monomer.atoms(force=True)]
        dssp.eval(monomer.paths["self"])
