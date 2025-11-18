import os
import sys

sys.path.append('..')

from src.bioiain.biopython import downloadPDB
import src.bioiain as bi

import numpy as np


bi.log("start", "test.py")

file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
                          file_path="./pdb_list.txt", file_format="pdb",
                          overwrite=False)

bi.log(1, "File folder:", file_folder)


structure = bi.imports.recover("5JJM")

if structure is None:
    structure = bi.imports.loadPDB(os.path.join(file_folder, "5JJM.pdb"))

bi.log("header", structure)

structure.init_all()

crystal = structure.get_crystals()

crystal.set_params(
    min_monomer_length=50,
    oligomer_levels=[2, 4],
    min_contacts=10,
    contact_threshold=15,
)

crystal.process(force="force" in sys.argv)



from src.bioiain.visualisation import pymol

script = pymol.PymolScript(folder=".", name="test")
script.load(crystal.paths["original"], "original", to="pdb")
script.cell()
script.symmetries()
script.group()
script.write_script()






bi.log("end", "DONE")
exit()


