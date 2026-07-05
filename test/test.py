import os, json, sys

import torch

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")
log("title", "test.py")



from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness_base import *

FOLDSEEK = os.environ.get("FOLDSEEK_PATH", "foldseek")


entity = CompactStructure.from_file("1M2Z.cif", export_folder="trash")
print(entity)


fs = FoldseekDB("db_" + entity.code(), [entity.path(source=True)], folder=entity.folder(), foldseek_command=FOLDSEEK)
for tensor, saprot_name, tensor_path, entry in fs.saprot_embeddings(save_folder=entity.folder()):
    code, chain, model = fs.parse_name(entry["name"])

    log("header", f"Code {code} chain {chain} model {model}")

    chains = entity.chains(chain, by_complex=True, model=model)
    for en in range(1280):
        entity.set_misc(f"{saprot_name}_{en}", None, force=False)
    for chain_entity in chains:
        print(chain_entity)
        residues = chain_entity.residues()
        print(len(residues), tensor.shape)

        assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}"


        chain_entity.img3D(property="b", plot=False, embedding=tensor)
