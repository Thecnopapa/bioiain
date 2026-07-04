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
set_seed()

IMG_SIZE = 16
INFERENCE = "-i" in sys.argv
TRAIN = "-t" in sys.argv

dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="pisa", oligo=0)
dataset.add_list("./data/cath-dataset-nonredundant-S20.multimeric.list", oligo=1)
dataset.shuffle()
dataset = StructureDataset.from_list("./data/consensus.list", name="consensus", oligo=1)

DATASET_NAME = f"{dataset.name}_3D_{IMG_SIZE}"

FORCE = False

embeddings = EmbeddingDataset(DATASET_NAME)
print(dataset)
if (not FORCE) or "--rebuild" in sys.argv:
    embeddings.load()

if (len(embeddings) == 0 or FORCE):
    for n, (entity, entry) in enumerate(dataset.entities(entity_class=CompactStructure, return_entries=True)):
        entity.compactness()
        for chain in entity.chains():
            img = chain.img3D(property="compactness", plot=False, size=16)

            exit()

            tensor = torch.tensor(img).reshape(-1, *img.shape)

            embedding = compactness3Dembedding.from_tensor(tensor,name=chain.full_id(), img_size=IMG_SIZE).save()
            key = embeddings.add(embedding, key=chain.full_id())
            embeddings.add_label(key, entry.oligo)

            #print(embeddings)
            #print(embeddings[len(embeddings)-1])
            #print(dataset)
            embeddings.save(temp=True)

    embeddings.save(temp=False)

if TRAIN:
    log("start", "TRAINING")
    model = Compactness3Dmk1(name=dataset.name, in_shape=[1, IMG_SIZE])
    model.mount()
    print(repr(model))

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
