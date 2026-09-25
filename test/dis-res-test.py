### src.bioiain import #################################################################################################
import sys
sys.path.append('..')
from src.bioiain import *
from src.bioiain.utilities import *
########################################################################################################################




dataset = StructureDataset.from_list(["6NWL"])
print(dataset)

for entity in dataset.entities():
    print(entity)
    print(entity.path())
    for res in entity.residues():
        print(res)