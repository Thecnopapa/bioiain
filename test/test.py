import os, json, sys


sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")


dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="monomers")

print(dataset)
exit()






from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness import *

entity = CompactStructure.from_file(os.path.join(".", "3sg0.cif"), export_folder="trash", force=False)
print(entity)


entity._calculate_compactness(session="--session" in sys.argv, plot="--plot" in sys.argv)







exit()




