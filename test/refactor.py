import os, json, sys
sys.path.append('..')


from src.bioiain.base import *

entity = BIEntity.from_file("./1M2Z.cif")
print(entity)
struc = entity.structure()
print(struc)



print(struc.chains())
for chain in struc.chains():

    print(chain.export())
    print(json.dumps(chain.data, indent=4))
    print(json.dumps(chain.paths, indent=4))

#print(struc.atoms())
#print(struc.residues())