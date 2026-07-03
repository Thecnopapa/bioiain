import os, json, sys

import torch

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")


from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness import *
from compactness_models import *

from src.bioiain.machine import *


FORCE = ("--force" in sys.argv) or ("-f" in sys.argv)
print(f"FORCE={FORCE}")



if "monomers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="monomers")

elif "multimers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.multimeric.list", name="multimers")

elif "pisa" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.multimeric.list", name="pisa")
    dataset.add_list("./data/cath-dataset-nonredundant-S20.multimeric.list")

elif "receptors" in sys.argv:
    dataset = StructureDataset.from_list("./data/receptors.list", name="receptors")

elif "lbds" in sys.argv:
    dataset = StructureDataset.from_list("./data/lbds.list", name="lbds")

elif "consensus" in sys.argv:
    dataset = StructureDataset.from_list("./data/consensus.list", name="consensus")

else:
    dataset = StructureDataset.from_list(["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"], name="aleph")


#dataset = StructureDataset.from_list("./data/consensus.list")
print(dataset)
#dataset.load(entity_class=CompactStructure, check_existing=not FORCE)
#dataset.export()



DATASET_NAME = f"{dataset.name}_foldseek"

log("start", "EMBEDDINGS")

embeddings = EmbeddingDataset(DATASET_NAME)
if not FORCE:
    embeddings.load()

if len(embeddings) == 0 or FORCE:
    n_missmatches = 0
    total = 0
    fs_cmd = os.environ.get("FOLDSEEK_PATH", "foldseek")
    fs = FoldseekDB(dataset.name, dataset, foldseek_command=fs_cmd)
    label_folder = os.path.join(SUBDIR_NAME, "labels", "compactness")
    os.makedirs(label_folder, exist_ok=True)
    for data, tensor in fs.match_dataset(dataset, atoms=False, entity_class=CompactStructure, force=FORCE):
        entity = data["entity"]
        chain = data["chain"]
        t = tensor[0]
        saprot_name = tensor[1]


        #print(t.shape)

        #print("DATA", data)
        #print("TENSOR", t.shape)
        #print("ENTITY:", entity)
        #print("CHAIN:", chain)

        #print(len(t), len(chain.sequence()))
        try:
            assert t.shape[-2] == len(chain.sequence()), f"Sequence len ({len(chain.sequence())}) and token len ({t.shape[-2]}) missmatch."
        except AssertionError as e:
            log("warning", e)
            n_missmatches += 1
            continue

        name = f"{chain.code()}_{chain.id()}"
        label_path = os.path.join(label_folder, f"{chain.code()}_{chain.id()}.compactness.label.csv")
        log(2, label_path)

        if FORCE or not os.path.exists(label_path):
            log(1, "Regenerating compactness label...")

            entity.compactness(with_symmetry=True, force=FORCE)
            chain = entity.chains(chain.complex(), by_complex=False)[0]
            #print("CHAIN 2:", chain)

            label_header="#"
            label_value=" "
            #print(len(chain.residues()))

            assert t.shape[-2] == len(chain.sequence()), f"Sequence len ({len(chain.sequence())}) and token len ({t.shape[-2]}) missmatch."

            for res in chain.residues():

                c = res.ca.get_misc("compactness")
                label_header += f"{res.resnum:8d} "
                if type(c) is float:
                   c = f"{c:7.3f}"
                label_value += f"{c}, "

            with open(label_path, "w") as lf:
                lf.write(label_header+"\n")
                lf.write(label_value[:-2]+"\n")
        else:
            log(1, "Compactness label already generated")
        log(2, label_path, saprot_name)

        embedding = SaProtEmbedding.from_tensor(t, name=name, subfolder=saprot_name).save()
        embeddings.add(embedding, label_path=label_path)
        embeddings.save(temp=True)

        total += 1
        print()



    log("header", f"Chains processed ok: {total-n_missmatches}, missmatches:{n_missmatches}")
    embeddings.save()

log("end", "EMBEDDINGS")

log("start", "TRAINING")

log("header", embeddings)


model = CompactnessMLPmk1(name=DATASET_NAME, in_shape=[1280], hidden_dims=[])
model.mount()
total_params = sum(p.numel() for p in model.submodels["default"].parameters())
log(1, "Number of parameters in the model:", total_params)
#print(model)

epochs = 100
for epoch in range(epochs):
    log("start", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)

    max_n = len(embeddings)
    for n, item in enumerate(embeddings):
        if n % 100 == 0:
            log(1, f"{n:6d}/{max_n:6d}", end = " ")
        if item.l is None:
            continue
        out = model.forward(item.t)
        #out = torch.clamp(out,0, 10)
        loss = model.loss(out, item)
        if n % 100 == 0:
            loss_str = f"{model.running_loss['default']/model.running_loss['total']:7.3f}"
            print(f"loss: {colour('yellow', loss_str)} \tlast: loss={loss:7.3f} out={out.item():7.3f} l={item.l:<7.3f}", end="\r")
    print()
    model.save(temp=True)
    model.add_epoch()
    log("end", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
model.save()



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
