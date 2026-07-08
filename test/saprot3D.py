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
if "--size" in sys.argv:
    IMG_SIZE = int(sys.argv[sys.argv.index("--size") + 1])
log(1, f"IMG_SIZE={IMG_SIZE}")

print(dataset)

MODEL_CLASS = Saprot3Dto1


def generate_3DSaprot_embeddings(dataset, img_size=16, foldseek_command=None, force=False, rebuild=False, force_labels=False, force_embeddings=False):
    from src.bioiain.machine.datasets import EmbeddingDataset
    if force:
        rebuild = True

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
            tracemalloc_top()
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
                assert len(chain) <= 1, f"Multiple chains detected {(code,ch,model)}: {chain}"
                chain = chain[0]
                log(1, chain)

                if not embedding_done:
                    residues = chain.residues(need_backbone=False)
                    log(1, "Loading tensor...")
                    try:
                        tensor = torch.load(tensor_path)
                    except Exception as e:
                        os.remove(tensor_path)
                        log("Error", "Error reading tensor:", tensor_path, e)
                        continue
                    #print(tensor.shape)
                    assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}"
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

            except (StructureLoadException, NotImplementedError) as e:
                dataset.add_to_blacklist(dataset.get(code).get("path"), e)
            except AssertionError as e:
                try:
                    log("warning", f'\n{entry["aa_seq"]}\n{chain.sequence()}')
                except:
                    pass
                dataset.add_to_blacklist(dataset.get(code).get("path"), e)

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


IN_SHAPE = embeddings.get(0).t.shape

if TRAIN:
    log("start", "TRAINING")
    log(1, "in_shape:", IN_SHAPE)

    model = MODEL_CLASS(name=dataset.name, in_shape=IN_SHAPE)
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
            out_i = model.forward(tensor)
            #print("OUT:", out_i.shape, out_c.shape)
            #print("LABEL:", label.shape, label.dtype)
            label = model.compress(label)
            #print("COMPRESSED LABEL:", label.shape)

            if label_oligo is not None:
                out_c = model.classify(out_i, reference=label)
                loss = model.raw_loss(out_i, label, out_c, label_oligo)
            else:
                loss = model.loss(out_i, label)
            #print("LOSS:", loss)

            if n % 1 == 0:
                loss_str = f"{model.running_loss['default'] / model.running_loss['total']:7.3f}"
                print(f"loss: {colour('yellow', loss_str)} \tlast --> loss={loss:7.3f} out={out_c.item():7.3f} l={item.l}",
                      end="\r")

            if "--preview" in sys.argv:
                from src.bioiain.visualisation import voxels3d
                out3d = out_i.detach().cpu().numpy()[0]
                count3d = (out3d != 0) & (out3d != 0) & (out3d != 0)
                voxels3d(out3d, count3d, show_plot=True, title=f"({item.name}) out={out_c.detach().item():5.3f} l={item.l}", shrink=True)
                exit()

        print()
        model.save(temp=True)
        model.add_epoch()
        log("end", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
    model.save()

    log("end", "TRAINING")


if INFERENCE:
    log("start", "INFERENCE")

    with torch.no_grad():
        try:
            filepath = sys.argv[sys.argv.index("--file") + 1]
            log(1, f"File path: {filepath}")
        except:
            raise


        try:
            model_path = sys.argv[sys.argv.index("--model") + 1]
        except:
            try:
                model_path = sys.argv[sys.argv.index("--md") + 1]
            except:
                model_path = None
        log(1, f"Model path: {model_path}")

        model = MODEL_CLASS(name=dataset.name, in_shape=IN_SHAPE, inference=True)
        log(1, f"Model:", model)
        model.load(model_path)


        log(1, "Loading entity...")
        entity = CompactStructure.from_file(filepath, export_folder="inference")
        entity.compactness(with_symmetry=True)

        entity.export()
        log(2, "Entity loaded:", entity)

        inference_name = f"inference_{MODEL_CLASS.__name__}_{datetime.datetime.now().strftime('_%y-%m-%d_%H-%M-%S')}"
        inference_folder = os.path.join(entity.folder(), "inference", inference_name)
        os.makedirs(inference_folder, exist_ok=False)



        fs = FoldseekDB("db_"+entity.code(), [entity.path(source=True)], folder=entity.folder())
        print(fs)

        for n, (tensor_path, entry, saprot_model) in enumerate(fs.saprot_embeddings(return_tensor=False)):
            code, ch, m = fs.parse_name(entry["name"])
            name = f"{code}_{ch}_{m}"

            try:
                log(1, entity)
                chain = entity.chains(ch, by_complex=True, model=m)
                print([c.complex() for c in chain])
                assert len(chain) <= 1, f"Multiple chains detected {(code,ch,model)}: {chain}"
                chain = chain[0]
                log(1, chain)

                residues = chain.residues(need_backbone=False)
                log(1, "Loading tensor...")
                try:
                    tensor = torch.load(tensor_path)
                except Exception as e:
                    os.remove(tensor_path)
                    log("Error", "Error reading tensor:", tensor_path, e)
                    continue
                #print(tensor.shape)
                assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}"
                log(1, "Generating 3D embedding...")
                tensor3D = chain.img3D(property=None, plot=False, size=IMG_SIZE, embedding=tensor, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, tensor3D.shape)


                log(1, "Loading compactness...")
                entity.compactness()
                log(1, "Generating 3D label...")
                label3D = chain.img3D(property="compactness", plot=False, size=IMG_SIZE, as_embedding=True, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, label3D.shape)

            except:
                raise

            in_tensor = tensor3D.to(DEVICE)
            print("IN_TENSOR", in_tensor.shape)
            real_label_x = label3D.to(DEVICE)
            print("REAL_LABEL", real_label_x.shape)
            compressed_label_x = model.compress(real_label_x)
            print("COMPRESSED_LABEL", compressed_label_x.shape)

            out_x = model(in_tensor)
            print("OUT", out_x.shape)

            import matplotlib.pyplot as plt
            from src.bioiain.visualisation import voxels3d, show
            fig = plt.figure()
            n_figs = 6

            label_ax = fig.add_subplot(1, n_figs, 1, projection="3d")
            real_label = real_label_x.detach().cpu().numpy()[0]
            label_count = (real_label > 0.1) & (real_label > 0.1) & (real_label > 0.1)
            label_count = label_count.astype(np.int64)
            #print(label_count)
            voxels3d(real_label, count_grid=label_count, ax = label_ax, shrink = True, title="Real label")

            compressed_ax = fig.add_subplot(1, n_figs, 2, projection="3d")
            compressed_label = compressed_label_x.detach().cpu().numpy()[0]
            voxels3d(compressed_label, ax = compressed_ax, shrink = True, title="Compressed")


            out_ax = fig.add_subplot(1, n_figs, 3, projection="3d")
            out = out_x.detach().cpu().numpy()[0]
            print(out)
            voxels3d(out, ax = out_ax, shrink = True, title="Out")

            scaled_ax = fig.add_subplot(1, n_figs, 4, projection="3d")

            max_val = compressed_label.reshape(compressed_label.shape[-1] ** 3).max()
            min_val = compressed_label.reshape(compressed_label.shape[-1] ** 3).min()

            max_out = out.reshape(out.shape[-1] ** 3).max()
            min_out = out.reshape(out.shape[-1] ** 3).min()

            out_range = abs(max_out - min_out)
            val_range = abs(max_val - min_val)
            
            scaled_out = np.absolute((out + min_out) / out_range * val_range)
            
            voxels3d(scaled_out, ax = scaled_ax, shrink = True, title="Scaled")


            diff_ax = fig.add_subplot(1, n_figs, 5, projection="3d")
            diff = np.absolute(scaled_out- compressed_label)
            voxels3d(diff, ax = diff_ax, shrink = True, title="Diff")

            model_diff_ax = fig.add_subplot(1, n_figs, 6, projection="3d")
            pred, model_diff = model.classify(out_x, reference=compressed_label_x, return_diff=True)
            pred = pred.detach().cpu().numpy().item()
            model_diff = model_diff.detach().cpu().numpy()[0]
            voxels3d(model_diff, ax = model_diff_ax, shrink = True, title=f"Model diff ({pred:.2f})")


            show()
            exit()


