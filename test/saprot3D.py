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


def generate_3DSaprot_embeddings(dataset, img_size=16, foldseek_command=None, force=False, rebuild=False, force_labels=False, force_embeddings=False):
    from src.bioiain.machine.datasets import EmbeddingDataset
    if force:
        rebuild = True

    if foldseek_command is None:
        foldseek_command = os.environ.get("FOLDSEEK_PATH", "foldseek")
    fs = FoldseekDB(dataset.name, dataset, foldseek_command=foldseek_command, force=force, dry=True)

    embeddings_name = f"{fs.saprot_model}_3D_size_{img_size}_{dataset.name}"
    labels_name = f"compactness_3D_size_{img_size}_{dataset.name}"

    embeddings = EmbeddingDataset(embeddings_name)
    labels = EmbeddingDataset(labels_name)

    if not force:
        if not force_embeddings:
            embeddings.load(load_temp=True)
        if not force_labels:
            labels.load(load_temp=True)

    if embeddings.incomplete() or labels.incomplete() or rebuild:
        fs.run()

        for n, (tensor_path, entry, saprot_model) in enumerate(fs.saprot_embeddings(return_tensor=False)):
            log(1, f"N={n}")
            code, ch, model = fs.parse_name(entry["name"])
            name = f"{code}_{ch}_{model}"
            embedding_done = False
            label_done = False
            if name in embeddings.embeddings.keys() and not force_embeddings:
                if os.path.exists(embeddings.embeddings[name]["embedding_path"]):
                    log(1, f"Embedding ({name}) already generated")
                    embedding_done = True
            if name in labels.embeddings.keys() and not force_labels:
                l_path = labels.embeddings[name]["embedding_path"]
                if os.path.exists(l_path):
                    log(1, f"Label ({name}) already generated")
                    label_done = True

            if (embedding_done and label_done) and not force:
                continue
            try:
                entity = CompactStructure.from_file(dataset.get(code).get("path"), verbose=False)
               
                log(1, entity)
                chain = entity.chains(ch, by_complex=True, model=model)
                print([c.complex() for c in chain])
                assert len(chain) == 1, f"Multiple chains detected {(code,ch,model)}: {chain}"
                chain = chain[0]
                log(1, chain)

                if not embedding_done:
                    residues = chain.residues(need_backbone=False)
                    log(1, "Loading tensor...")
                    tensor = torch.load(tensor_path)
                    #print(tensor.shape)
                    assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}\n{entry["aa_seq"]}\n{chain.sequence()}"
                    log(1, "Generating 3D embedding...")
                    tensor3D = chain.img3D(property=None, plot=False, size=IMG_SIZE, embedding=tensor, mode="mean", residue_kwargs={"need_backbone":False})
                    log(2, tensor3D.shape)
                    embedding = SaProt3DEmbedding.from_tensor(tensor3D,name=name, img_size=IMG_SIZE, saprot_model=saprot_model).save()
                    embeddings.add(embedding)
                    embeddings.save(temp=True)

                if not label_done:
                    log(1, "Loading compactness...")
                    entity.compactness()
                    log(1, "Generating 3D label...")
                    label3D = chain.img3D(property="compactness", plot=False, size=IMG_SIZE, as_embedding=True, mode="mean", residue_kwargs={"need_backbone":False})
                    log(2, label3D.shape)
                    label_embedding = Compactness3DEembedding.from_tensor(label3D, name=name, img_size=IMG_SIZE).save()
                    labels.add(label_embedding)
                    labels.save(temp=True)

            except StructureLoadException as e:
                dataset.add_to_blacklist(dataset.get(code).get("path"), e)
            except AssertionError:
                raise
    embeddings.save(temp=False)
    labels.save(temp=False)
    return embeddings, labels


embeddings, labels = generate_3DSaprot_embeddings(dataset, img_size=IMG_SIZE, force=FORCE, rebuild=REBUILD, force_labels=LABELS)

if REBUILD or FORCE or LABELS:
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
    in_shape = embeddings.get(0).t.shape
    log(1, "in_shape:", in_shape)

    model = Saprot3Dto1(name=dataset.name, in_shape=in_shape)
    print(model)
    model.mount()
    print(repr(model))


    EPOCHS = 100
    for epoch in range(EPOCHS):
        log("start", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
        max_n = len(embeddings)
        for n in range(len(embeddings)):
            #print(n)
            embeddings.use_label("oligo")
            item = embeddings.get(n, label=False, label_key="oligo")
            litem = labels.get(n, label=False)
            assert item.name == litem.name
            tensor = item.t.to(torch.float32).to(DEVICE)
            label = litem.t.to(torch.float32).to(DEVICE)
            label_oligo = item.l
            #print(tensor, label_oligo)

            #print(entity)
            if n % 1 == 0:
                log(1, f"{n:6d}/{max_n:6d}", end=" ")


            #print("\nIN:", tensor.shape, tensor.dtype)
            out_i, out_c = model.forward(tensor)
            #print("OUT:", out_i.shape, out_c.shape)
            #print("LABEL:", label.shape, label.dtype)
            label = model.compress(label)
            #print("COMPRESSED LABEL:", label.shape)
            #print(out_i is None, label is None, out_c is None, label_oligo is None)
            if label_oligo is not None:
                loss = model.raw_loss(out_i, label, out_c, label_oligo)
            else:
                loss = model.loss(out_i, label)
            #print("LOSS:", loss)

            if n % 1 == 0:
                loss_str = f"{model.running_loss['default'] / model.running_loss['total']:7.3f}"
                print(f"loss: {colour('yellow', loss_str)} \tlast --> loss={loss:7.3f} out={out_c.item():7.3f} l={item.l}",
                      end="\r")

            continue # Remove to plot example output
            from src.bioiain.visualisation import voxels3d
            out3d = out_i.detach().cpu().numpy()[0]
            count3d = (out3d != 0) & (out3d != 0) & (out3d != 0)
            voxels3d(out3d, count3d, show_plot=True, title=f"({item.name}) out={out_c.detach().item():5.3f} l={item.l}", shrink=True)


        print()
        model.save(temp=True)
        model.add_epoch()
        log("end", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
    model.save()

    log("end", "TRAINING")


if INFERENCE:
    pass
