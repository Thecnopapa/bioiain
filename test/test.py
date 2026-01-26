import os, sys
sys.path.append('..')



from src.bioiain.biopython import downloadPDB, Structure
from src.bioiain.symmetries import Crystal
import src.bioiain as bi

import numpy as np

force = "force" in sys.argv

bi.log("start", "test.py")

#file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
#                          file_path="./pdb_list.txt", file_format="pdb",
#                          overwrite=False)
file_folder = downloadPDB("/home/iain/projects/vib-ai/internship/data", "receptors",
                          file_path="/home/iain/projects/vib-ai/internship/data/receptors.txt", file_format="cif",
                          overwrite=False)

bi.log(1, "File folder:", file_folder)

for n, file in enumerate(os.listdir(file_folder)):
    if not file.endswith(".cif"):
        continue
    code = file[:4]
    #structure = bi.biopython.recover(code)
    bi.log("title", code)


    if force:
        structure = bi.biopython.loadPDB(os.path.join(file_folder, f"{code}.cif"))
        structure.init_all()
    else:
        structure = Structure.recover(code, data_path=f"exports/{code}/{code}", load_structure=True)

    bi.log("header", structure)

    if force:
        crystals = structure.init_crystal()


    crystals = structure.data["crystals"]
    print(crystals, structure.paths["crystal_folder"])
    crystal = Crystal.recover("cryst", data_path=os.path.join(structure.paths["crystal_folder"], crystals[0]), load_structure=True)

    crystal.set_crystal_params(
        min_monomer_length=50,
        min_contacts=10,
        contact_threshold=15,
    )


    if force:
        if crystal.process(force=force) is None:

            exit()
    monomers = crystal.data["monomers"]
    print("monomers")
    print(monomers)
    from src.bioiain.symmetries.interactions import *
    for monomer in monomers:
        interactions_per_monomer(monomer, crystal.paths["monomer_folder"])




    exit()
    continue

    # crystal.get_oligomers(
    #     oligomer_levels=[2],
    # )


    # from src.bioiain.symmetries import Oligomer
    #
    # print(crystal.paths["oligo_folder"])
    # for file in sorted(os.listdir(crystal.paths["oligo_folder"]))   :
    #     print(file)
    #     if file.endswith(".data.json"):
    #         oligo = Oligomer.recover(id="recovered",data_path=os.path.join(crystal.paths["oligo_folder"], file))
    #         print(oligo)
    #
    #
    # from src.bioiain.visualisation import pymol
    #
    # script = pymol.PymolScript(folder=".", name="test")
    # script.load(crystal.paths["original"], "original", to="pdb")
    # script.cell()
    # script.symmetries()
    # script.group()
    # script.write_script()






bi.log("end", "DONE")
exit()


