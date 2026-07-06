import os, json, sys

import torch

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")
log("title", "test.py")



from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness_base import *
from compactness_models import *


FOLDSEEK = os.environ.get("FOLDSEEK_PATH", "foldseek")


entity = CompactStructure.from_file("1M2Z.cif", export_folder="trash")
print(entity)



IMG_SIZE = 16
fs = FoldseekDB("db_" + entity.code(), [entity.path(source=True)], folder=entity.folder(), foldseek_command=FOLDSEEK)
for tensor, saprot_name, tensor_path, entry in fs.saprot_embeddings(save_folder=entity.folder()):
    code, chain, model = fs.parse_name(entry["name"])

    log("header", f"Code {code} chain {chain} model {model}")

    chains = entity.chains(chain, by_complex=True, model=model)
    for chain_entity in chains:
        print(chain_entity)
        residues = chain_entity.residues()
        print(len(residues), tensor.shape)

        assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}"


        tensor3D = chain_entity.img3D(property=None, plot=False, size=IMG_SIZE, embedding=tensor, mode="mean")
        print(tensor3D.shape)

        embedding = SaProt3DEmbedding.from_tensor(tensor,name=chain_entity.full_id(), img_size=IMG_SIZE).save()
