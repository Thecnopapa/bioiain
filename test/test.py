import os, json, sys, time, random
sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.base import *

log("start", "test.py")
log("title", "test.py")


from src.bioiain.symmetry import *
from src.bioiain.symmetry.niggli import *


cell = UnitCell(3., 5.196, 2., 103.55, 109.28, 134.53)
cell = UnitCell(3., 5.196, 2., 103.92, 109.47, 134.88)

print(cell)
print(cell.is_buerger(), cell.is_niggli())
niggli = cell_to_niggli(cell, verbose=False, eps=0.01)
print(niggli)
print(niggli.is_buerger(), niggli.is_niggli())
exit()

structures = StructureDataset.from_list(["1OII", "1VZ4"], name="symmetry_test",)
print(structures)


for entity in structures.entities():
    log("header", entity)
    log(1, entity.space_group())
    cell = entity.unit_cell()
    log(1, cell)
    niggli = cell_to_niggli(cell)
    log(1, niggli)






log("end")
