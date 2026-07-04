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
entity.compactness()
fs = FoldseekDB("db_" + entity.code(), [entity.path(source=True)], folder=entity.folder(), foldseek_command=FOLDSEEK)
for tensor, saprot_name, tensor_path, entry in fs.saprot_embeddings(save_folder=entity.folder()):
    code, chain, model = fs.parse_name(entry["name"])

    log("header", f"Code {code} chain {chain}")

    chain_entity = entity.chains(chain, use_complex=True, model=model)[0]
    print(chain)
    residues = chain_entity.residues()

    chain_entity.img3D(property=["b", "compactness"], plot=True)
