import os, json, sys
sys.path.append('..')


from src.bioiain.base import *
from src.bioiain.base.mmcif import *

entity = BIEntity.from_file("./3HHB.cif")
print(entity)
struc = entity.structure()
print(struc)



print(struc.chains())
for chain in struc.chains():

    print(chain.export())
    print(json.dumps(chain.data, indent=4))
    print(json.dumps(chain.paths, indent=4))

    write_dict(chain.data["info"], "bi_info", file_path="testdict", name=chain.name())

#print(struc.atoms())
#print(struc.residues())


for n, atom in enumerate(struc.atoms(hetatm=True, disordered=True)):
    if False:
        continue
    print(atom.pdb_string())

struc.export(as_pdb=True)
