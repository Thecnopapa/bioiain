

import os
import sys

sys.path.append('..')
from src.bioiain.utilities import log
from src.bioiain.biopython.imports import *
from src.bioiain.symmetries.crystal import *
from src.bioiain.biopython.SASA import SASA





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
        log("start", "SASA")
        print(monomer)
        print()
        
        surfece_res_ids = monomer.get_surface_residues()
        print(surfece_res_ids)
        monomer.export(include_misc=True)
        log("end", "SASA")


