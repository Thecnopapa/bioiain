import os, json, sys


sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")


from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness import *




dataset = StructureDataset.from_list("./data/consensus.list")
dataset.load(entity_class=CompactStructure)
print(dataset)

#fs_cmd = "/cri4/iain/bin/foldseek/bin/foldseek"
fs_cmd="foldseek"

fs = FoldseekDB(dataset.name, dataset, foldseek_command=fs_cmd)

for data, tensor in fs.match_dataset(dataset, atoms=False, entity_class=CompactStructure):
    print("DATA", data)
    print("TENSOR", tensor[1])
    print("ENTITY:", data["entity"])


print("DONE")
exit()




entity = BIEntity.from_file(os.path.join(".", "3sg0.cif"), export_folder="trash", force=True)
print(entity)

if "v2" in sys.argv:

    entity.fragment(aleph_mode="ALEPH2", in_place=True)
elif "original" in  sys.argv:
    entity.fragment(in_place=True)
else:
    entity.fragment(aleph_mode="debug", in_place=True)
exit()

entity._calculate_compactness(session="--session" in sys.argv, plot="--plot" in sys.argv)




exit()



