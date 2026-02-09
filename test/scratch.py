

import os
import sys

sys.path.append('..')
from src.bioiain.utilities import log
from src.bioiain.biopython.imports import *
from src.bioiain.symmetries.crystal import *
from src.bioiain.biopython.SASA import SASA
from src.bioiain.visualisation.pymol import *





file_folder = downloadPDB("./data",
                          list_name="mari-weird",
                          file_format="cif",
                          pdb_list=["1M2Z"],
                          overwrite=False)

sasa = SASA()
for  file in os.listdir(file_folder):
    if "1M2Z" not in file:
        continue
    monomers = get_monomers(file, file_folder)



    for monomer in monomers:
        monomer.show_exposed_residues(force=True)


