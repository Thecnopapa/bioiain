import os
import sys

sys.path.append('..')

from src.bioiain.biopython import downloadPDB
import src.bioiain as bi

import numpy as np


bi.log("start", "test.py")

#file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
#                          file_path="./pdb_list.txt", file_format="pdb",
#                          overwrite=False)
file_folder = downloadPDB("/home/iain/projects/vib-ai/internship/data", "receptors",
                          file_path="/home/iain/projects/vib-ai/internship/data/receptors.txt", file_format="cif",
                          overwrite=False)

bi.log(1, "File folder:", file_folder)

for file in os.listdir(file_folder):
    if not file.endswith(".cif"):
        continue
    code = file[:4]
    #structure = bi.biopython.recover(code)
    bi.log("title", code)
    structure = None
    if structure is None:
        structure = bi.biopython.loadPDB(os.path.join(file_folder, f"{code}.cif"))

    bi.log("header", structure)

    structure.init_all()

    crystal = structure.get_crystals()

    crystal.set_params(
        min_monomer_length=50,
        oligomer_levels=[2],
        min_contacts=10,
        contact_threshold=15,
    )

    if crystal.process(force="force" in sys.argv) is None:
        continue


    from src.bioiain.symmetries import Oligomer

    print(crystal.paths["oligo_folder"])
    for file in sorted(os.listdir(crystal.paths["oligo_folder"]))   :
        print(file)
        if file.endswith(".data.json"):
            oligo = Oligomer.recover(id="recovered",data_path=os.path.join(crystal.paths["oligo_folder"], file))
            print(oligo)


    from src.bioiain.visualisation import pymol

    script = pymol.PymolScript(folder=".", name="test")
    script.load(crystal.paths["original"], "original", to="pdb")
    script.cell()
    script.symmetries()
    script.group()
    script.write_script()






bi.log("end", "DONE")
exit()


