import os, json, sys

import torch

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "saprot3D.py")
log("title", "saprot3D.py")



from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness_base import *
from compactness_models import *
from src.bioiain.machine import *
from data import dataset
set_seed()



IMG_SIZE = 16

print(dataset)


def generate_3DSaprot_embeddings(dataset, img_size=16, foldseek_command=None, force=False, rebuild=False):
    from src.bioiain.machine.datasets import EmbeddingDataset
    if force:
        rebuild = True

    if foldseek_command is None:
        foldseek_command = os.environ.get("FOLDSEEK_PATH", "foldseek")
    fs = FoldseekDB(dataset.name, dataset, foldseek_command=foldseek_command, force=force, dry=True)

    embeddings_name = f"{fs.saprot_model}_3D_size_{img_size}_{dataset.name}"
    embeddings = EmbeddingDataset(embeddings_name)
    if not force:
        embeddings.load(load_temp=True)

    if embeddings.incomplete() or rebuild:
        fs.run()

        for n, (tensor_path, entry, saprot_model) in enumerate(fs.saprot_embeddings(return_tensor=False)):
            log(1, f"N={n}")
            code, ch, model = fs.parse_name(entry["name"])
            name = f"{code}_{ch}_{model}"
            if name in embeddings.embeddings.keys():
                log(1, f"Embedding ({name}) already generated")
                continue
            try:
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
            except StructureLoadException as e:
                dataset.add_to_blacklist(dataset.get(code).get("path"), e)
            except AssertionError:
                raise
    embeddings.save(temp=False)
    return embeddings


embeddings = generate_3DSaprot_embeddings(dataset, img_size=IMG_SIZE, force=FORCE, rebuild=REBUILD)

if REBUILD or FORCE:
    log("header","Configuring oligomer labels")
    for n, k in enumerate(embeddings.embeddings.keys()):
        log(2, f"{n+1}/{len(embeddings)}", end="\r")
        code = k.split("_")[0]
        entry = dataset.get(code)
        label = entry.get("oligo", None)
        embeddings.add_label(k, label, "oligo")
    embeddings.use_label("oligo")
    embeddings.save()
    print()
    log(1, "Oligomer labels ready")

if TRAIN:
    log("start", "TRAINING")

    in_shape = embeddings[0].t.shape
    log(1, "in_shape:", in_shape)

    model = Saprot3Dto1(name=dataset.name, in_shape=in_shape)
    print(model)
    exit()
    model.mount()
    print(repr(model))
    exit()

    EPOCHS = 100
    for epoch in range(EPOCHS):
        log("start", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
        max_n = len(dataset)
        for n, item in enumerate(embeddings):

            #print(entity)
            if n % 100 == 0:
                log(1, f"{n:6d}/{max_n:6d}", end=" ")

            out = model.forward(item.t.reshape(-1, *item.t.shape).to(DEVICE))
            loss = model.loss(out, item)
            if n % 100 == 0:
                loss_str = f"{model.running_loss['default'] / model.running_loss['total']:7.3f}"
                print(f"loss: {colour('yellow', loss_str)} \tlast: loss={loss:7.3f} out={out.item():7.3f} l={entry.oligo}",
                      end="\r")
        print()
        model.save(temp=True)
        model.add_epoch()
        log("end", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
    model.save()

    log("end", "TRAINING")


if INFERENCE:
    pass
