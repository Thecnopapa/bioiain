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
print(t)
print(t.id)

t.init_crystal()
t.export()
t.pass_down()

model = t.get_list()[0]
print(model)
print(model.export())
print(model.get_full_id(), model.data["info"]["name"])


crystal = bi.symmetries.Crystal.cast(model.copy())
crystal.set_params(
    min_monomer_length=50,
    oligomer_levels=[2],
)


crystal.process()
crystal.export()

print(crystal, model)



from src.bioiain.visualisation import pymol


script = pymol.PymolScript(folder=".", name="test")
script.load(t.paths["original"], "original", to="pdb")
script.cell()
script.symmetries()
script.group()
script.write_script()







exit()

script._bioiain = "sys\nsys.path.append('..')\nimport src.bioiain"


script.print("pdb")
script.disable("(all)")
script.add("bi.log", "'header'", "\"I'm a log\"", is_cmd = False)






t.export("./exports", data=True)


script.execute()









