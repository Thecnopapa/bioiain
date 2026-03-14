import os, json, sys
sys.path.append('..')


from src.bioiain.base import *


struc = BIStructure().from_file("./1M2Z.cif")

print(struc.atoms())
print(struc.residues())