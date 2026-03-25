import os, json, sys

sys.path.append('..')

from src.bioiain.utilities import *
from src.bioiain.utilities.logging import *

log("start", "ligands.py")

tracemalloc_start()

from src.bioiain.base import *



try:
	entity = FragmentedStructure.recover_from_id("1M2Z", export_folder="./a")
except:
	entity = BIEntity.from_file("./1M2Z.cif", export_folder="./a")

log("header", "ligands")

for ligand in entity.ligands():
	log(1,ligand)
	for a in ligand.atoms:
		#log(2, a)
		pass


entity = entity.fragment()
print(entity.cvmatrix())

entity.export()



log("end", "DONE")
