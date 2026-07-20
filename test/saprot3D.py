import os, json, sys, time, threading

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
from src.bioiain.utilities.parallel import N_THREADS
from data import dataset
set_seed()


ALLOW_EXPORTS= not ("--no-exports" in sys.argv)
log(1, f"ALLOW_EXPORTS={ALLOW_EXPORTS}")

IMG_SIZE = 8
if "--size" in sys.argv:
    IMG_SIZE = int(sys.argv[sys.argv.index("--size") + 1])
log(1, f"IMG_SIZE={IMG_SIZE}")

log(1, dataset)


MODEL_CLASS = Saprot3Dto1

WORK_AS_IS = "--as-is" in sys.argv


def generate_3DSaprot_embeddings(dataset, 
                                 img_size=16, 
                                 foldseek_command=None, 
                                 force=False, rebuild=False, 
                                 force_labels=False, 
                                 force_embeddings=False, 
                                 as_is=False, 
                                 allow_exports=True, 
                                 n_threads=1):
    from src.bioiain.machine.datasets import EmbeddingDataset
    if force:
        rebuild = True
        allow_exports = False

    fs = FoldseekDB(dataset.name, dataset, foldseek_command=foldseek_command, force=force, dry=True)

    embeddings_name = f"{fs.saprot_model}_3D_size_{img_size}_{dataset.name}"
    rel_labels_name = f"compactness_SYM_3D_size_{img_size}_{dataset.name}"
    abs_labels_name = f"compactness_ABS_3D_size_{img_size}_{dataset.name}"

    embeddings = EmbeddingDataset(embeddings_name)
    rel_labels = EmbeddingDataset(rel_labels_name)
    abs_labels = EmbeddingDataset(abs_labels_name)

    if not force:
        if not force_embeddings or as_is:
            embeddings.load(load_temp=True)
        if not force_labels or as_is:
            rel_labels.load(load_temp=True)
            abs_labels.load(load_temp=True)

    log(1, "Loaded datasets:")
    log(2, embeddings)
    log(2, rel_labels)
    log(2, abs_labels)

    def generate_embedding_set(n, 
                               name, 
                               code, 
                               ch, 
                               model, 
                               tensor_path,
                               entry, 
                               saprot_model, 
                               img_size, 
                               embedding_done,
                               label_done,
                               ):
        try:
            entity = FragmentedStructure.from_file(dataset.get(code).get("path"), verbose=False, check_existing=allow_exports)
           
            log(1, entity)
            chain = entity.chains(ch, by_complex=True, model=model)
            try:
                assert len(chain) <= 1, f"Multiple chains detected {(code,ch,model)}: {chain}"
            except:
                print([c.complex() for c in chain])
                raise MultipleChainsDetected(f"Multiple chains detected {(code,ch,model)}: {chain}")

            try:
                chain = chain[0]
            except:
                print(entity.chains(by_complex=True, model=model))
                raise NoChainsDetected(f"No chains detected {(code,ch,model)}: {chain} {print(entity.chains(by_complex=True, model=model))}")

            log(1, chain)

            if not embedding_done:
                residues = chain.residues(need_backbone=False)
                log(1, "Loading tensor...")
                try:
                    tensor = torch.load(tensor_path)
                except Exception as e:
                    os.remove(tensor_path)
                    log("Error", "Error reading tensor:", tensor_path, e)
                    return
                #print(tensor.shape)
                try:
                    assert tensor.shape[-2] == len(residues), f"{tensor.shape[-2]} / {len(residues)}"
                except:
                    log("warning", f'\n{entry["aa_seq"]}\n{chain.sequence()}')
                    raise SequenceMissmatchException(f"{tensor.shape[-2]} / {len(residues)}")
                log(1, "Generating 3D embedding...")
                tensor3D = chain.img3D(property=None, plot=False, size=img_size, embedding=tensor, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, tensor3D.shape)
                embedding = SaProt3DEmbedding.from_tensor(tensor3D,name=name, img_size=img_size, saprot_model=saprot_model).save()
                embeddings.add(embedding)
                embeddings.save(temp=True)
                log(2, embeddings)

            if not label_done:
                log(1, "Loading Relative compactness...")
                entity.compactness(with_symmetry=True)
                log(1, "Generating Relative 3D label...")
                rel_label3D = chain.img3D(property="rel_compactness", plot=False, size=img_size, as_embedding=True, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, rel_label3D.shape)

                log(1, "Loading Absolute compactness...")
                chain.compactness(with_symmetry=False, export=False)
                log(1, "Generating Absolute 3D label...")
                abs_label3D = chain.img3D(property="abs_compactness", plot=False, size=IMG_SIZE, as_embedding=True, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, abs_label3D.shape)

                log(1, "Adding labels to datasets...")

                rel_label_embedding = Compactness3DEembedding.from_tensor(rel_label3D, name=name, img_size=img_size, relative=True).save()
                rel_labels.add(rel_label_embedding)
                rel_labels.save(temp=True)
                log(2, rel_labels)

                abs_label_embedding = Compactness3DEembedding.from_tensor(abs_label3D, name=name, img_size=IMG_SIZE, relative=False).save()
                abs_labels.add(abs_label_embedding)
                abs_labels.save(temp=True)
                log(2, abs_labels)

        except (StructureLoadException, NotImplementedError, MultipleChainsDetected, NoChainsDetected, SequenceMissmatchException, ALEPHError) as e:
            dataset.add_to_blacklist(dataset.get(code).get("path"), e)
        except AssertionError as e:
            try:
                log("warning", f'\n{entry["aa_seq"]}\n{chain.sequence()}')
            except:
                pass
            dataset.add_to_blacklist(dataset.get(code).get("path"), e)

    if (embeddings.incomplete() or rel_labels.incomplete() or abs_labels.incomplete() or rebuild) and not as_is:
        fs.run()

        threads = []
        n_total = len(fs)
        for n, (tensor_path, entry, saprot_model) in enumerate(fs.saprot_embeddings(return_tensor=False)):
            tracemalloc_top()
            log(1, f"N={n+1}/{n_total}")
            code, ch, model = fs.parse_name(entry["name"])
            name = f"{code}_{ch}_{model}"
            embedding_done = False
            label_done = False
            if name in embeddings.embeddings.keys() and not force_embeddings:
                if os.path.exists(embeddings.embeddings[name]["embedding_path"]):
                    log(1, f"Embedding ({name}) already generated")
                    embedding_done = True
            if name in rel_labels.embeddings.keys() and name in abs_labels.embeddings.keys() and not force_labels:
                rel_l_path = rel_labels.embeddings[name]["embedding_path"]
                abs_l_path = abs_labels.embeddings[name]["embedding_path"]
                if os.path.exists(rel_l_path) and os.path.exists(abs_l_path):
                    log(1, f"Labels ({name}) already generated")
                    label_done = True

            if (embedding_done and label_done) and not force:
                continue
            if N_THREADS > 1:
                while threading.active_count() > n_threads:
                    print(f"Waiting for available thread... ({name}) running:{threading.active_count()-1}")
                    time.sleep(1)
                thread = threading.Thread(target=generate_embedding_set,
                                          name = name,
                                          kwargs = dict(n=n,
                                                       name=name, 
                                                       code=code, 
                                                       ch=ch, 
                                                       model=model, 
                                                       tensor_path=tensor_path, 
                                                       entry=entry, 
                                                       saprot_model=saprot_model, 
                                                       img_size=img_size, 
                                                       embedding_done=embedding_done,
                                                       label_done=label_done, 
                                                       )
                                          )
                threads.append(t)
                thread.start()
            else:
                generate_embedding_set(n=n,
                                       name=name, 
                                       code=code, 
                                       ch=ch, 
                                       model=model, 
                                       tensor_path=tensor_path, 
                                       entry=entry, 
                                       saprot_model=saprot_model, 
                                       img_size=img_size, 
                                       embedding_done=embedding_done,
                                       label_done=label_done, 
                                       )
        
        while threading.active_count() > 1:
            print(f"waiting for Threads to finish")
            time.sleep(1)
        print("Joining threads...")
        for t in threads:
            t.join()
        print("All threads joined")


        embeddings.save(temp=False)
        rel_labels.save(temp=False)
        abs_labels.save(temp=False)
    return embeddings, rel_labels, abs_labels


embeddings, rel_labels, abs_labels = generate_3DSaprot_embeddings(dataset, img_size=IMG_SIZE, force=FORCE, allow_exports=ALLOW_EXPORTS, rebuild=REBUILD, force_embeddings=EMBEDDINGS, force_labels=LABELS, as_is=WORK_AS_IS, n_threads=N_THREADS)
log(1, "EMBEDDINGS", embeddings)
log(1, "REL LABELS:", rel_labels)
log(1, "ABS_LABELS:", abs_labels)
print(embeddings.keys()[-10:])
print(rel_labels.keys()[-10:])
print(abs_labels.keys()[-10:])
rel_labels.sort_as(embeddings)
abs_labels.sort_as(embeddings)
print(embeddings.keys()[-10:])
print(rel_labels.keys()[-10:])
print(abs_labels.keys()[-10:])
log(1, "EMBEDDINGS", embeddings)
log(1, "REL LABELS:", rel_labels)
log(1, "ABS_LABELS:", abs_labels)

if REBUILD or FORCE or LABELS:
    log("header","Configuring oligomer labels")
    for n, k in enumerate(embeddings.embeddings.list()):
        log(2, f"{n+1}/{len(embeddings)}", end="\r")
        code = k.split("_")[0]
        label = None
        if not "--precalculated" in sys.argv:
            with open(f"./data/{dataset.name}.monomeric.list") as mf:
                for l in mf:
                    if code in l:
                        label = 0
                        break
            if label is None:
                with open(f"./data/{dataset.name}.multimeric.list") as mf:
                    for l in mf:
                        if code in l:
                            label = 1
                            break
        else:
            entry = dataset.get(code)
            label = entry.get("oligo", None)
        embeddings.add_label(k, label, "oligo")
    embeddings.use_label("oligo")

    embeddings.save(temp=WORK_AS_IS)
    print()
    log(1, "Oligomer labels ready")
    print()


MODEL_NAME = dataset.name

DIFFERENCES = ("--diffs" in sys.argv) or ("--differences" in sys.argv) or ("--diff" in sys.argv)

if DIFFERENCES:
    MODEL_NAME+"_diffs"

IN_SHAPE = embeddings.get(0).t.shape
log(1, f"IN_SHAPE={IN_SHAPE}")

log(1, "EMBEDDINGS", embeddings)
log(1, "REL LABELS:", rel_labels)
log(1, "ABS_LABELS:", abs_labels)
PRINT_EVERY=100
SAVE_TEMP_MODELS= False
FINETUNE = "--finetune" in sys.argv
FAST = "--fast" in sys.argv
if len(embeddings) <= 1000 or DEVICE == "cpu":
    PRINT_EVERY=1
    if not FAST:
        SAVE_TEMP_MODELS = True

if TRAIN:
    log("start", "TRAINING")

    model = MODEL_CLASS(name=MODEL_NAME, in_shape=IN_SHAPE)
    log(1, model)
    model.mount()
    print(repr(model))


    EPOCHS = 100
    for epoch in range(EPOCHS):
        log("start", f"EPOCH: {epoch}", print_timer=True, reset_timer=False)
        max_n = len(embeddings)
        n = 0
        for n in range(len(embeddings)):
            #print(n)
            embeddings.use_label("oligo")
            try:
                item = embeddings.get(n, label=True, label_key="oligo")
            except DeletedIndex:
                continue
            key = item.key
            rel_item = rel_labels.get(rel_labels.get_indexes(key), label=False)
            abs_item = abs_labels.get(abs_labels.get_indexes(key), label=False)

            assert item.name == rel_item.name and item.name == abs_item.name, f"{item.name} == {rel_item.name} == {abs_item.name}"
            tensor = item.t.to(torch.float32).to(DEVICE)
            rel_label = rel_item.t.to(torch.float32).to(DEVICE)
            abs_label = abs_item.t.to(torch.float32).to(DEVICE)
            label_oligo = torch.Tensor([item.l]) if item.l is not None else None
            #print(n, item.name, label_oligo, item.l)

            #print(entity)
            if n % PRINT_EVERY == 0:
                log(1, f"{n+1:6d}/{max_n:6d}", end=" ")

            #print("\nIN:", tensor)
            #print("\nIN SHAPE:", tensor.shape, tensor.dtype)
            out_i = model.forward(tensor)
            #print("OUT:", out_i)
            #print("OUT SHAPE:", out_i.shape)
            #print(rel_item.t)
            #print("REL LABEL:", rel_label)
            #print("REL LABEL SHAPE:", rel_label.shape, rel_label.dtype)
            rel_label_x = model.compress(rel_label)
            abs_label_x = model.compress(abs_label)

            if DIFFERENCES:
                label = torch.subtract(rel_label_x, abs_label_x)
            else:
                label = rel_label_x
            #print("COMPRESSED REL LABEL:", rel_label_x)
            #print("COMPRESSED REL LABEL:", rel_label_x.shape)

            if label_oligo is not None and FINETUNE:
                out_c = model.classify(out_i, reference=label)
                loss = model.raw_loss(out_i, label, out_c, label_oligo)
                out_c_text = f"{out_c.item():7.3f}"
            else:
                out_c_text="None"
                loss = model.loss(out_i, label)
            #print("LOSS:", loss)

            if n % PRINT_EVERY == 0:
                loss_str = f"{model.running_loss['default'] / model.running_loss['total']:15.3f}"
                print(f"loss: {colour('yellow', loss_str)}\tlast --> loss={loss.item():15.3f} out={out_c_text} l={item.l}",
                      end="\r")

            if "--preview" in sys.argv:
                from src.bioiain.visualisation import voxels3d
                out3d = out_i.detach().cpu().numpy()[0]
                count3d = (out3d != 0) & (out3d != 0) & (out3d != 0)
                voxels3d(out3d, count3d, show_plot=True, title=f"({item.name}) out={out_c.detach().item():5.3f} l={item.l}", shrink=True)
                exit()

        print()
        if epoch % 10 == 0 and epoch != 0:
            if SAVE_TEMP_MODELS:
                model.save(temp=True)

        model.add_epoch(add_histograms=(epoch % 10 == 0) and not FAST)
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

        model = MODEL_CLASS(name=MODEL_NAME, in_shape=IN_SHAPE, inference=True)
        log(1, f"Model:", model)
        model.load(model_path)


        log(1, "Loading entity...")
        entity = FragmentedStructure.from_file(filepath, export_folder="inference", check_existing=False)
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
                log(1, entity, name)
                chain = entity.chains(ch, by_complex=True, model=m)
                print([c.complex() for c in entity.chains()])
                print(chain, len(chain))
                assert len(chain) <= 1, f"Multiple chains detected {(code,ch,m)}: {chain}"
                assert len(chain) > 0, f"No chains detected {(code,ch,m)}: {entity.chains(by_complex=True)}"
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


                log(1, "Loading Relative compactness...")
                entity.compactness(with_symmetry=True, force=True)
                log(1, "Generating relative 3D label...")
                rel_label3D = chain.img3D(property="rel_compactness", shrink=True, size=IMG_SIZE, as_embedding=True, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, rel_label3D.shape)

                log(1, "Loading Absolute compactness...")
                chain.compactness(with_symmetry=False, export=False)
                log(1, "Generating absolute 3D label...")
                abs_label3D = chain.img3D(property="abs_compactness", shrink=True, size=IMG_SIZE, as_embedding=True, mode="mean", residue_kwargs={"need_backbone":False})
                log(2, abs_label3D.shape)

            except Exception as e:
                log("warning", e)
                continue
                raise

            in_tensor = tensor3D.to(DEVICE)
            rel_label = rel_label3D.to(DEVICE)
            rel_label_x = model.compress(rel_label)

            abs_label = abs_label3D.to(DEVICE)
            abs_label_x = model.compress(abs_label)

            out_x = model(in_tensor)
            #print("OUT", out_x.shape)

            import matplotlib.pyplot as plt
            from src.bioiain.visualisation import voxels3d, show, close
            n_figs = 6
            n_rows = 3
            fig = plt.figure(figsize=(n_figs*5, n_rows*5))
            fig.suptitle(f"{chain}")



            # Real Labels ##############################################################################################
            rel_label_ax = fig.add_subplot(n_rows, n_figs, 1, projection="3d")
            rel_label_detached = rel_label.detach().cpu().numpy()[0]
            rel_label_count = (rel_label_detached > 0.01) & (rel_label_detached > 0.01) & (rel_label_detached > 0.01)
            rel_label_count = rel_label_count.astype(np.int64)
            #print(rel_label_detached)
            voxels3d(rel_label_detached, count_grid=rel_label_count, ax = rel_label_ax, shrink = True, title="Rel label")

            abs_label_ax = fig.add_subplot(n_rows, n_figs, n_figs+1, projection="3d")
            abs_label_detached = abs_label.detach().cpu().numpy()[0]
            abs_label_count = (abs_label_detached > 0.01) & (abs_label_detached > 0.01) & (abs_label_detached > 0.01)
            abs_label_count = abs_label_count.astype(np.int64)
            #print(abs_label_detached)
            voxels3d(abs_label_detached, count_grid=abs_label_count, ax = abs_label_ax, shrink = True, title="Abs label")


            diff_label_ax = fig.add_subplot(n_rows, n_figs, n_figs*2+1, projection="3d")

            diff_label = np.subtract(abs_label.detach().cpu().numpy()[0],rel_label.detach().cpu().numpy()[0])

            diff_label_count = (abs(diff_label) > 0.001) & (abs(diff_label) > 0.001) & (abs(diff_label) > 0.001)
            diff_label_count = diff_label_count.astype(np.int64)
            #print(diff_label)
            #print(diff_label_count)
            #print(diff_label.shape)
            voxels3d(diff_label, count_grid=diff_label_count, ax = diff_label_ax, shrink = True, title="Diff label")


            # Compressed labels ########################################################################################
            compressed_rel_ax = fig.add_subplot(n_rows, n_figs, 2, projection="3d")
            rel_label_x_detached = rel_label_x.detach().cpu().numpy()[0]
            voxels3d(rel_label_x_detached, ax = compressed_rel_ax, shrink = True, title="Rel label X")

            compressed_abs_ax = fig.add_subplot(n_rows, n_figs, n_figs+2, projection="3d")
            abs_label_x_detached = abs_label_x.detach().cpu().numpy()[0]
            voxels3d(abs_label_x_detached, ax = compressed_abs_ax, shrink = True, title="Abs label X")

            compressed_diff_ax = fig.add_subplot(n_rows, n_figs, n_figs*2+2, projection="3d")
            diff_label_x = abs(rel_label_x.detach().cpu().numpy()[0] - abs_label_x.detach().cpu().numpy()[0])
            voxels3d(diff_label_x, ax = compressed_diff_ax, shrink = True, title="Diff label X")



            # Model output #############################################################################################
            row = 1
            if DIFFERENCES:
                row=2
            out_ax = fig.add_subplot(n_rows, n_figs, n_figs*row+3, projection="3d")
            out = out_x.detach().cpu().numpy()[0]
            #print(out)
            voxels3d(out, ax = out_ax, shrink = True, title="Raw output")

            scaled_ax = fig.add_subplot(n_rows, n_figs, n_figs*row+4, projection="3d")

            max_val = rel_label_x_detached.reshape(rel_label_x_detached.shape[-1] ** 3).max()
            min_val = rel_label_x_detached.reshape(rel_label_x_detached.shape[-1] ** 3).min()

            max_out = out.reshape(out.shape[-1] ** 3).max()
            min_out = out.reshape(out.shape[-1] ** 3).min()

            #print(max_out, min_out)
            #print(max_val, min_val)
            out_range = abs(max_out - min_out)
            val_range = abs(max_val - min_val)

            #print(out_range)
            #print(val_range)
            
            #scaled_out = np.absolute((out + min_out) / out_range * val_range)
            scaled_out = ((out + min_out) / out_range) * val_range

            voxels3d(scaled_out, ax = scaled_ax, shrink = True, title="Scaled output")

            # Real vs output differences ################################################################################
            rel_diff_ax = fig.add_subplot(n_rows, n_figs, 5, projection="3d")

            rel_label_x_detached = rel_label_x.detach().cpu().numpy()[0]

            max_rel_vals = np.maximum.reduce([scaled_out, rel_label_x_detached])
            min_rel_vals = np.minimum.reduce([scaled_out, rel_label_x_detached])

            rel_diff = np.absolute(max_rel_vals - min_rel_vals)
            voxels3d(rel_diff, ax = rel_diff_ax, shrink = True, title="Diff out/rel")

            abs_diff_ax = fig.add_subplot(n_rows, n_figs, n_figs+5, projection="3d")

            abs_label_x_detached = abs_label_x.detach().cpu().numpy()[0]

            max_abs_vals = np.maximum.reduce([scaled_out, abs_label_x_detached])
            min_abs_vals = np.minimum.reduce([scaled_out, abs_label_x_detached])

            abs_diff = np.absolute(max_abs_vals - min_abs_vals)
            voxels3d(abs_diff, ax = abs_diff_ax, shrink = True, title="Diff out/abs")


            model_rel_diff_ax = fig.add_subplot(n_rows, n_figs, 6, projection="3d")
            rel_pred, model_rel_diff = model.classify(out_x, reference=rel_label_x, return_diff=True)
            rel_pred = rel_pred.detach().cpu().numpy().item()
            model_rel_diff_detached = model_rel_diff.detach().cpu().numpy()[0]
            voxels3d(model_rel_diff_detached, ax = model_rel_diff_ax, shrink = True, title=f"Model diff out/rel (p:{rel_pred:.2f})")

            model_abs_diff_ax = fig.add_subplot(n_rows, n_figs, n_figs+6, projection="3d")
            abs_pred, model_abs_diff = model.classify(out_x, reference=abs_label_x, return_diff=True)
            abs_pred = abs_pred.detach().cpu().numpy().item()
            model_abs_diff_detached = model_abs_diff.detach().cpu().numpy()[0]
            voxels3d(model_abs_diff_detached, ax = model_abs_diff_ax, shrink = True, title=f"Model diff out/abs (p:{abs_pred:.2f})")



            # DIFf diffs
            diff_diff_ax = fig.add_subplot(n_rows, n_figs, n_figs*2+5, projection="3d")
            max_diff_diff_vals = np.maximum.reduce([rel_diff, abs_diff])
            min_diff_diff_vals = np.minimum.reduce([rel_diff, abs_diff])
            diff_diff = np.absolute(max_diff_diff_vals - min_diff_diff_vals)

            diff_diff_count = (diff_diff >= 0.001) & (diff_diff > 0.001) & (diff_diff > 0.001)
            diff_diff_count = diff_diff_count.astype(np.int64)


            voxels3d(diff_diff, count_grid=diff_diff_count, ax = diff_diff_ax, shrink = True, title="Diff diff")






            img_path = os.path.join(inference_folder, name+".png")
            plt.savefig(img_path, dpi=300)
            log(1, f"Figure saved to:", img_path)
            show()
            close(fig)
            #exit()

