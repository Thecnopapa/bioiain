import os, json, sys

import torch

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "compactness.py")


from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness_base import *
from compactness_models import *

from src.bioiain.machine import *


set_seed()



FORCE = ("--force" in sys.argv) or ("-f" in sys.argv)
print(f"FORCE={FORCE}")



if "monomers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="monomers")

elif "multimers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.multimeric.list", name="multimers")

elif "pisa" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="pisa")
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


ENTITY_CLASS = CompactStructure

DATASET_NAME = f"{dataset.name}_foldseek"
log("title", DATASET_NAME)

IN_SHAPE=[1280]
FOLDSEEK = os.environ.get("FOLDSEEK_PATH", "foldseek")


log("start", "EMBEDDINGS")

embeddings = EmbeddingDataset(DATASET_NAME)
if (not FORCE) or "--rebuild" in sys.argv:
    embeddings.load()

if len(embeddings) == 0 or FORCE:
    n_missmatches = 0
    total = 0
    fs = FoldseekDB(dataset.name, dataset, foldseek_command=FOLDSEEK)
    label_folder = os.path.join(SUBDIR_NAME, "labels", "compactness")
    os.makedirs(label_folder, exist_ok=True)
    for n, (data, tensor) in enumerate(fs.match_dataset(dataset, atoms=False, entity_class=ENTITY_CLASS, force=FORCE)):
        log("header", n)
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
            #log("warning", e)
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


MODEL_CLASS = CompactnessMLPmk1
log("header", f"MODEL_CLASS={MODEL_CLASS}")


if "-t" in sys.argv:
    log("start", "TRAINING")

    epochs = 100


    log(1, embeddings)
    log(1, f"EPOCHS={epochs}")


    model = MODEL_CLASS(name=DATASET_NAME, in_shape=IN_SHAPE)
    model.mount()
    total_params = sum(p.numel() for p in model.submodels["default"].parameters())
    log(1, "Number of parameters in the model:", model.n_params(human=True))
    print(repr(model))
    #print(model)

    for epoch in range(epochs):
        log("start", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)

        max_n = len(embeddings)
        for n, item in enumerate(embeddings):
            if n % 100 == 0:
                log(1, f"{n:6d}/{max_n:6d}", end = " ")
            if item.l is None:
                continue
            out = model.forward(item.t.to(DEVICE))
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

    log("end", "TRAINING")



if "-i" in sys.argv:
    log("start", "INFERENCE")

    with torch.no_grad():
        filepath = sys.argv[sys.argv.index("--file") + 1]
        log(1, f"File path: {filepath}")
        try:
            model_path = sys.argv[sys.argv.index("--model") + 1]
        except:
            model_path = None
        log(1, f"Model path: {model_path}")


        log(1, "Loading entity...")
        entity = ENTITY_CLASS.from_file(filepath, export_folder="inference")
        entity.export()
        log(2, "Entity loaded:", entity)

        inference_name = f"inference_{MODEL_CLASS.__name__}_{datetime.datetime.now().strftime('_%y-%m-%d_%H-%M-%S')}"
        inference_folder = os.path.join(entity.folder(), "inference", inference_name)
        os.makedirs(inference_folder, exist_ok=False)

        fs = FoldseekDB(entity.code(), [entity.path(source=True)], folder=entity.folder(), foldseek_command=FOLDSEEK)
        print(fs)



        model = MODEL_CLASS(name=DATASET_NAME, in_shape=IN_SHAPE, inference=True)
        log(1, f"Model:", model)
        model.load(model_path)





    log("end", "INFERENCE")


print("DONE")
exit()