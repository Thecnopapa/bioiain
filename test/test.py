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


from data import dataset
IMG_SIZE = 16

def generate_3DSaprot_embeddings(dataset, img_size=16, foldseek_command=None, force=False, rebuild=False):
    from src.bioiain.machine.datasets import EmbeddingDataset
    if force:
        rebuild = True

    if foldseek_command is None:
        foldseek_command = os.environ.get("FOLDSEEK_PATH", "foldseek")
    fs = FoldseekDB(dataset.name, dataset, foldseek_command=foldseek_command, force=force, dry=True)

    embeddings_name = f"{fs.saprot_model}_3D_size_{img_size}_{dataset.name}"
    embeddings = EmbeddingDataset(embeddings_name)
    if not rebuild:
        embeddings.load(load_temp=True)

    if embeddings.incomplete():
        fs.run()

        for n, (tensor_path, entry, saprot_model) in enumerate(fs.saprot_embeddings(return_tensor=False)):
            log(1, f"N={n}")
            code, ch, model = fs.parse_name(entry["name"])
            name = f"{code}_{ch}_{model}"
            if name in embeddings.embeddings.keys():
                log(1, f"Embedding ({name}) already generated")
                continue
            entity = BIEntity.from_file(dataset.get(code).get("path"), verbose=False)
            log(1, entity)
            chain = entity.chains(ch, by_complex=True, model=model)
            assert len(chain) == 1, f"Multiple chains detected {(code,ch,model)}: {chain}"
            chain = chain[0]
            log(1, chain)
            residues = chain.residues(need_backbone=False)
            log(1, "Loading tensor...")
            tensor = torch.load(tensor_path)
            #print(tensor.shape)
            assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}\n{entry["aa_seq"]}\n{chain.sequence()}"
            log(1, "Generating 3D embedding...")
            tensor3D = chain.img3D(property=None, plot=False, size=IMG_SIZE, embedding=tensor, mode="mean", residue_kwargs={"need_backbone":False})
            #print(tensor3D.shape)
            embedding = SaProt3DEmbedding.from_tensor(tensor,name=name, img_size=IMG_SIZE, saprot_model=saprot_model).save()
            embeddings.add(embedding)
            embeddings.save(temp=True)
    embeddings.save(temp=False)
    return embeddings

embeddings = generate_3DSaprot_embeddings(dataset, img_size=IMG_SIZE)
print(embeddings)

exit()





fs = FoldseekDB("db_" + entity.code(), [entity.path(source=True)], folder=entity.folder(), foldseek_command=FOLDSEEK)
for tensor, saprot_name, tensor_path, entry in fs.saprot_embeddings():
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
