import os
import sys



sys.path.append('..')
import Bio.PDB as bp



from src.bioiain.biopython import downloadPDB
from src.bioiain.biopython.base import BiopythonOverlayClass
import src.bioiain as bi



file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
                          file_path="./pdb_list.txt", file_format="pdb",
                          overwrite=False)

bi.log("header", file_folder)


t = bi.imports.loadPDB(os.path.join(file_folder,os.listdir(file_folder)[0]))
#print(t)
print(t.id)

t.init_crystal()
t.export("./exports", data=True, structure=False)

model = t.get_list()[0]
print(model)

crystal = bi.symmetries.Crystal.cast(model.copy())
crystal.set_params(
    data = t.data,
    min_monomer_length=100,
    oligomer_levels=[2],
)


crystal.process()
crystal.export_data("./exports", "crystal")

print(crystal, model)












exit()
from src.bioiain.visualisation import pymol


script = pymol.PymolScript()
script._bioiain = "sys\nsys.path.append('..')\nimport src.bioiain"

script.load(t.paths["original"], "original", to="pdb")
script.print("pdb")
script.disable("(all)")
script.add("bi.log", "'header'", "\"I'm a log\"", is_cmd = False)






t.export("./exports", data=True)

script.write_script("./exports")
script.execute()









