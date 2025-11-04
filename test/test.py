import os
import sys
sys.path.append('.')
import Bio.PDB as bp



from src.bioiain.biopython import downloadPDB
from src.bioiain.biopython.base import BiopythonOverlayClass
import src.bioiain as bi



file_folder = downloadPDB("./test/data", "test_list", ["5JJM", "6nwl"],
                          file_path="test/pdb_list.txt", file_format="pdb",
                          overwrite=False)

bi.log("header", file_folder)


t = bi.imports.loadPDB(os.path.join(file_folder,os.listdir(file_folder)[0]))
print(t)
print(t.id)

t.init_crystal()

print(t.export("./test/exports", data=True))










