import os, json, sys, time, random
sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.base import *

log("start", "test.py")
log("title", "test.py")


from src.bioiain.symmetry import *

structures = StructureDataset.from_list(["1OII", "1VZ4"], name="symmetry_test",)
print(structures)


for entity in structures.entities():
    print(entity)
    print(entity.space_group())
    print(entity.unit_cell())






log("end")
