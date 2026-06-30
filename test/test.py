import os, json, sys


sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")


from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness import *




dataset = StructureDataset.from_list("./data/consensus.list")
print(dataset)

#fs_cmd = "/cri4/iain/bin/foldseek/bin/foldseek"
fs_cmd="foldseek"

fs = FoldseekDB(dataset.name, dataset, foldseek_command=fs_cmd)


fasta = fs.tokens_fasta()

for name, t_seq in fasta.parse().items():
    name = name.split(" ")[0]
    print(name.split("_"))
    if len(name.split("_")) == 1:
        code = name.split("_")[0]
        chain = "*"
    elif len(name.split("_")) == 2:
        code, chain = name.split("_")
    else:
        raise Exception(f"Unable to fetch name and code from {name}")
    print(code, chain)
    t_seq = t_seq[0]
    print(code, chain, len(t_seq))
    
    entry = dataset.get(code)
    print(entry)

exit()
print(fs.tokens_fasta().get_names())




for entity in dataset.entities():
    print(entity)
    for chain in entity.chains():
        print(chain.sequence())
    print()

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



