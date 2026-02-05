

import os
import sys

sys.path.append('..')

from src.bioiain.biopython.imports import *
from src.bioiain.symmetries.crystal import *


file_folder = downloadPDB("./data",
                          list_name="mari-weird",
                          file_format="cif",
                          pdb_list=["1QO1", "1QCR", "3M9C", "1FE1", "1M2Z"],
                          overwrite=False)

print(file_folder)
print(os.listdir(file_folder))


for  file in os.listdir(file_folder):
    monomers = get_monomers(file, file_folder)
    print(monomers)
    for monomer in monomers:
        [print(a) for a in monomer.atoms(force=True)]