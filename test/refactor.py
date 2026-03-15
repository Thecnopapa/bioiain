import os, json, sys
sys.path.append('..')


from src.bioiain.base import *
from src.bioiain.base.mmcif import *

entity = BIEntity.from_file("./3HHB.cif")
print(entity)
struc = entity.structure()
struc.export()
print(struc)
#print(struc.headers)

print(struc.params())
print(struc._get_operations())