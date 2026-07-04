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
from src.bioiain.machine import *
set_seed()

IMG_SIZE = 16

dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="pisa", oligo=0)
dataset.add_list("./data/cath-dataset-nonredundant-S20.multimeric.list", oligo=1)
#dataset.shuffle()



log("start", "TRAINING")
model = Compactness3Dmk1(name=dataset.name, in_shape=[1, IMG_SIZE])
model.mount()
print(repr(model))



EPOCHS = 100
for epoch in range(EPOCHS):
    log("start", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
    max_n = len(dataset)
    for n, (entity, entry) in enumerate(dataset.entities(entity_class=CompactStructure, return_entries=True)):
        #entity = BIEntity.from_file("1M2Z.cif", export_folder="trash")
        entity.compactness()
        entity.export()
        #print(entity)
        if n % 100 == 0:
            log(1, f"{n:6d}/{max_n:6d}", end=" ")

        for chain in entity.chains():
            #print(chain)


            img = chain.img3D(property="compactness", plot=False, size=16)

            np.set_printoptions(threshold=sys.maxsize)
            #print(img)
            #exit()

            tensor = torch.tensor(img).reshape(-1, *img.shape).to(DEVICE)
            #print(tensor)
            #print(tensor.shape)


            out = model.forward(tensor)
            #print(out)
            label = torch.tensor([[float(entry.oligo)]]).to(DEVICE)
            #print(label)
            loss = model.loss(out, label)
            #print(loss)
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
