import os, sys, json
sys.path.append('..')
from src.bioiain.biopython import downloadPDB
from src.bioiain import log
from src.bioiain.utilities.parallel import *
from src.bioiain.utilities.exceptions import *
import asyncio
import datetime


log("title", "test.py")
log("start", "SET UP")



BLACKLIST = ["1JJ2"]

#file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
#                          file_path="./pdb_list.txt", file_format="pdb",
#                          overwrite=False)
if "--no-download" in sys.argv:
    if "monomers" in sys.argv:
        file_folder = "./data/cath-monomeric"
        pdb_list = "cath"
    elif "cath" in sys.argv:
        file_folder = "./data/cath-nonredundant-S20"
        pdb_list = "cath"
    else: 
        file_folder = "./data/receptors"
        pdb_list = "cath"
else:
    if "monomers" in sys.argv:
        file_folder = downloadPDB("./data", "cath-monomeric",
                                               file_path="./data/cath-dataset-nonredundant-S20.monomeric.list",
                                               file_format="cif",
                                               overwrite=False)
        pdb_list = "monomers"

    elif "cath" in sys.argv:
        file_folder = downloadPDB("./data", "cath-nonredundant-S20",
                                               file_path="./data/cath-dataset-nonredundant-S20.list",
                                               file_format="cif",
                                               overwrite=False)
        pdb_list = "cath"
    else:
        file_folder = downloadPDB("./data", "receptors",
                                  file_path="./data/receptors.txt", file_format="cif",
                                  overwrite=False)
        pdb_list="rcps"

    log("header", "File folder:", file_folder)


from src.bioiain.symmetries.elements import Monomer
from src.bioiain.symmetries.crystal import get_monomers
from src.bioiain.machine.datasets import EmbeddingDataset
from src.bioiain.machine.embeddings import SaProtEmbedding, MissingProgram, FoldseekError
from src.bioiain.symmetries.interactions import InteractionProfile
from src.bioiain.utilities.sequences import MSA

import torch
torch.set_num_threads(avail_cpus)
log(1, f"Torch using {avail_cpus} threads")



#mlog = mem_log()

from src.bioiain.machine import DEVICE


FORCE = "--force" in sys.argv or "-f" in sys.argv

THRESHOLD=10
if "--threshold" in sys.argv:
    THRESHOLD = int(sys.argv[sys.argv.index("--threshold") + 1])

REBUILD = "--rebuild" in sys.argv
if FORCE:
    pass
ONLY = None
if "--only" in sys.argv:
    ONLY = sys.argv[sys.argv.index("--only") + 1].split(",")
DUAL = True
DUAL_CLASSES = True
DUAL_CLASSES_V2 = True
DUAL_CLASSES_V3 = False

if "--v1" in sys.argv:
    DUAL_CLASSES_V2 = False
elif "--v3" in sys.argv:
    DUAL_CLASSES_V2 = False
    DUAL_CLASSES_V3 = True

data_name = f"saprot_interactions_{pdb_list}_T{THRESHOLD}"
if DUAL:
    data_name += "_DUAL"
    if DUAL_CLASSES:
        data_name += "_CLASSES"
        if DUAL_CLASSES_V2:
            data_name += "_V2"
        elif DUAL_CLASSES_V3:
            data_name += "_V3"



if ONLY is not None:
    data_name += "_ONLY_" + "_".join(ONLY)

log("header", f"Data name: {data_name}")

dataset = EmbeddingDataset(name=data_name)
if not REBUILD:
    dataset.load()

from src.bioiain.machine import models


MODEL_NAME = "Golden"

if "--mk" in sys.argv: 
    MODEL_NAME = sys.argv[sys.argv.index("--mk") + 1]

else:
    if "mk3" in sys.argv:
        MODEL_CLASS = MLP_MK3
    elif "mk4" in sys.argv:
        MODEL_CLASS = DUAL_MLP_MK4
    elif "mk5" in sys.argv:
        MODEL_CLASS = DUAL_MLP_MK5
    elif "mk6" in sys.argv:
        MODEL_CLASS = DUAL_MLP_MK6
    elif "mk7" in sys.argv:
        MODEL_CLASS = DUAL_MLP_MK7
    elif "mk8" in sys.argv:
        MODEL_CLASS = DUAL_MLP_MK8
    elif "mk9" in sys.argv:
        MODEL_CLASS = DUAL_MLP_MK9

MODEL_CLASS = getattr(models, MODEL_NAME)

log(1, f"Model class: {MODEL_CLASS}")

LR = 0.0001
if "--lr" in sys.argv:
    LR = float(sys.argv[sys.argv.index("--lr") + 1])

log(1, f"Learning rate: {LR}")

MIX = False
if "--mix" in sys.argv:
    MIX = True

SKIP_ABS = False
if "--skip-abs" in sys.argv:
    SKIP_ABS = True

SHORT_EPOCHS = False
if "--short-e" in sys.argv:
    SHORT_EPOCHS = True

BATCH_SIZE = 32
if "--batch" in sys.argv:
    BATCH_SIZE = int(sys.argv[sys.argv.index("--batch") + 1])
log(1, "BATCH_SIZE:", BATCH_SIZE)

SEED = None
if "--seed" in sys.argv:
    SEED = int(sys.argv[sys.argv.index("--seed") + 1])
log(1, "SEED:", SEED)


#Comment

if "-l" in sys.argv or "-e" in sys.argv:

    if (not "absolute_calcuated" in dataset.data or REBUILD) and not SKIP_ABS:
        log("start", "ABSOLUTE LABELS")

        n_files = len(os.listdir(file_folder))
        for n, file in enumerate(sorted(os.listdir(file_folder))):

            if not file.endswith(".cif"):
                continue
            if any([b in file for b in BLACKLIST]):
                log("warning", f"File: {file} in blacklist!")
                continue
            if os.path.getsize(os.path.join(file_folder, file)) > 1.5 * 1024 * 1024:
                log("warning", f"File: {file} too large! ({os.path.getsize(os.path.join(file_folder, file)) * 1024 * 1024} MiB)")
                continue

            if ONLY is not None:
                if file[:4] not in ONLY:
                    continue
            log("title", f"ABSOLUTE {n}/{n_files}")
            try:
                mon_data = get_monomers(file, file_folder, only_ids=True, force=FORCE, contact_threshold=15)
            except StructureRecoverException:
                print("retrying")
                mon_data = get_monomers(file, file_folder, only_ids=True, force=True, contact_threshold=15)
                monomer = Monomer.recover(data_path=os.path.join(dataset.data["export_folder"], monomer_id.split("_")[0], "monomers", monomer_id))
            except Exception as e:
                log("warning", "Error processing crystal (skipped): {e}")
                continue
            if mon_data is None:
                log("Warning", f"{file} has no monomers")
                continue

            monomers, monomer_folder = mon_data

            if dataset.data.get("export_folder", None) is None:
                dataset.data["export_folder"] = "./exports"
            for monomer_id in monomers:
                try:
                    log("header", f"Calculating absolute interactions for: {monomer_id}")
                    if monomer_id in dataset and not FORCE and not REBUILD:
                        log(1, f"{monomer_id}: already in dataset")
                        continue

                    monomer = Monomer.recover(data_path=os.path.join(dataset.data["export_folder"], monomer_id.split("_")[0], "monomers", monomer_id))

                    if monomer is None:
                        log("Warning", f"{monomer_id} has no monomer")
                        continue
                    log(1, "Generating embeddings...")
                    embedding = SaProtEmbedding(entity=monomer, force=FORCE)
                    key = dataset.add(embedding=embedding, key=monomer.get_name(), fasta=True)
                    log(1, "Generating absolute labels...")
                    ints = InteractionProfile(monomer, threshold=THRESHOLD, force=FORCE)
                    label = ints.generate_labels(relative=False, force=FORCE)
                    dataset.add_label_from_string(label, key=key, var_name="abs_label")
                except MissingProgram as e:
                    raise e

                except FoldseekError as e:
                    log("warning", e)
                    #mon_data = get_monomers(file, file_folder, only_ids=True, force=True, contact_threshold=15)
                    continue
                except Exception as e:
                    log("Error", f"Exception occurred processing: {monomer_id}:\n", e)
                    raise e

                dataset.save()


        dataset.data["absolute_calcuated"] = True
        dataset.save()




    if (not "relative_calcuated" in dataset.data) or REBUILD or "--rel" in sys.argv:
        log("start", "RELATIVE LABELS")

        if DUAL_CLASSES:
            dataset.use_label("abs_label")
        dataset.map()


        datset_path = dataset.save()
        log("header", "DATASET:", dataset)
        msa = MSA(dataset.data["fasta_path"], dataset.data["name"], verbose=True)
        log("header", "MSA:", msa)

        n_mons = len(list(dataset.embeddings.keys()))
        for n, monomer_id in enumerate(list(dataset.embeddings.keys())):
            log("title", f"RELATIVE {n}/{n_mons}")
            try:
                log("header", f"Calculating relative interactions for: {monomer_id}")

                if "rel_label" in dataset.embeddings[monomer_id] and not (FORCE or REBUILD or "-rel" in sys.argv):
                    if dataset.embeddings[monomer_id]["rel_label"] is not None :
                        log(1, f"{monomer_id}: relative interactions already calculated")
                        continue

                monomer = Monomer.recover(data_path=os.path.join(dataset.data["export_folder"], monomer_id.split("_")[0], "monomers", monomer_id))

                if monomer is None:
                    log("Warning", f"{monomer_id} has no monomer")
                    exit()
                    continue
                log(1, "Generating relative labels...")
                ints = InteractionProfile(monomer, threshold=THRESHOLD, force=FORCE)

                rel_label, n_neighbours = ints.generate_labels(relative=True, force=FORCE, dataset=dataset, msa=msa, dual=DUAL, in_lab_var="abs_label", V2=DUAL_CLASSES_V2, V3=DUAL_CLASSES_V3)
                print("REL LAB:", len(rel_label), f"DUAL={DUAL}")

                if DUAL:
                    if DUAL_CLASSES:
                        if DUAL_CLASSES_V2:
                            dataset.add_label_from_list(rel_label, key=monomer_id, var_name="dual_discrete_2")
                        elif DUAL_CLASSES_V3:
                            dataset.add_label_from_list(rel_label, key=monomer_id, var_name="dual_discrete_3")

                        else:
                            dataset.add_label_from_list(rel_label, key=monomer_id, var_name="dual_class_label")

                        dataset.embeddings[monomer_id]["n_neighbours"] = n_neighbours

                    else:
                        dataset.add_label_from_list(rel_label, key=monomer_id, var_name="dual_label")
                else:
                    dataset.add_label_from_list(rel_label, key=monomer_id, var_name="rel_label")
                print(dataset)


                dataset.save()

            except SequenceMissmatchException:
                log("warning", "Sequence missmatch in:", monomer_id)
                dataset.remove(monomer_id)
                dataset.save()
                continue


        dataset.data["relative_calcuated"] = True
        dataset.save()



if "-t" in sys.argv:
    log("title", "TRAINING")
    log("start", "TRAINING")

    log("header", f"Dataset: {dataset}")
    if DUAL_CLASSES_V2:
        dataset.use_label("dual_discrete_2")
    elif DUAL_CLASSES_V3:
        dataset.use_label("dual_discrete_3")
    elif DUAL_CLASSES:
        dataset.use_label("dual_class_label")
    elif DUAL:
        dataset.use_label("dual_label")
    else:
        dataset.use_label("rel_label")
    dataset.split(seed=SEED)
    if DUAL and not DUAL_CLASSES:
        label_to_index = dataset.map(label_to_index={"contactability":0, "outer":1})
    else:
        label_to_index = dataset.map()
    log(2, "Label map:")
    print(json.dumps(label_to_index, indent=4))
    label_count = dataset.data["lab_count"]
    log(2, "Label count:")
    print(json.dumps(label_count, indent=4))

    epochs = 10
    if "--epochs" in sys.argv:
        epochs = int(sys.argv[sys.argv.index("--epochs") + 1])


    run_name = f"{data_name}_LR{str(LR).split(".")[-1]}_E{epochs}_MIX{int(MIX)}_B{BATCH_SIZE}"
    model = MODEL_CLASS(name=run_name, in_shape=(1280,), num_classes=len(label_to_index), lr=LR, weights=label_count, batch_size=BATCH_SIZE).to(DEVICE)

    model.add_map(dataset)
    run_name = model.data["name"]
    model.add_histogram("relative_neighbours", [int(e.get("n_neighbours", 0)) for e in dataset.embeddings.values() if not e.get("deleted", False)])
    model.add_text("hparams", json.dumps({
        "model_name": model.__class__.__name__,
        "dataset": data_name,
        "label": dataset.data["label_key"],
        "seed": SEED,
        "optimiser": model.optimisers["default"].__class__.__name__,
        "loss_fn": model.criterions["default"].__class__.__name__,
        "lr": LR,
        "batch_size": BATCH_SIZE,
        "mix": MIX,
        "target_epochs": epochs,
        "device": DEVICE,
        }, indent=4))



    log(1, f"Run name: {run_name}")
    total_params = sum(p.numel() for p in model.submodels["default"].parameters())
    log(1, "Number of parameters in the model:", total_params)



    log("header", f"Epochs: {epochs}")
    data_len = len(dataset)
    model.add_epoch()
    for epoch in range(epochs):
        epoch = epoch +1

        log("header", "EPOCH:", epoch)
        if MIX:
            dataset.split()

        dataset.train()


        if "no-dropout" in model.layers.keys() and epoch > 20:
            model.set_mode("no-dropout")
            log(1, "Disabling dropout layers...")
            print(repr(model))

        for optimiser in model.optimisers.values():
            log(1, "Learning rates")
            for group in optimiser.param_groups:
                log(2, group['lr'], type(group['lr']))

        if not is_cluster:
            epoch_start = datetime.datetime.now()
            eta = datetime.timedelta(0)
        print()
        try:
            for n, item in enumerate(dataset):
                if SHORT_EPOCHS:
                    if n > 100: break
                if not is_cluster:
                    print(f"\033]0;TRAIN {epoch}/{epochs} {(n/len(dataset))*100:3.0f}%\a", end="\r")
                #print("tensor:", item.t)
                #print("label:", item.l)
                #print(DEVICE)
                item.to(DEVICE)
                out = model(item.t)
                #print("out:", out)
                #print(type(item.l), item.l)

                loss = model.loss(out, item)
                #print(out.item(), item.l, loss)
                #exit()
                #model.step()

                if not is_cluster:
                    if n % (data_len // 100) == 0:
                        now = datetime.datetime.now()
                        delta = now - epoch_start
                        eta = datetime.timedelta(seconds=(delta.total_seconds() / (n+1)) * (data_len - n))
                    print(f"{n:6d}/{len(dataset):6d}: LOSS={model.running_loss['default']/model.running_loss['total']:6.4f} ({loss:6.4f}) <-- o:{torch.max(out, dim=0)[1].item()} t:{item.li} ETA: {eta}", end = "\r")

            model.save()
            model.test(dataset)
            model.send_run(host="iainvisa.com", key=os.environ.get("IAINVISA_FILE_KEY", None), epoch=epoch)
            model.add_epoch()

        except KeyboardInterrupt:
            log("warning", "!! Stopping model...")
            try:
                #print(os.environ)
                #print("KEY:", os.getenv("IAINVISA_FILE_KEY"))
                #model.send_run(host="iainvisa.com", key=os.environ.get("IAINVISA_FILE_KEY"))
                #print(model.send_run(host="127.0.0.1:5000", key=os.getenv("IAINVISA_FILE_KEY"), protocol="http").__dict__)
                print(repr(model))
                model.test(dataset, temp=True)
            except ModelNotFound as e:
                print(e)
            break





SCRIPT = None
TEMP = "--temp" in sys.argv
if "-p" in sys.argv:
    log("title", "PREDICTING...")
    from src.bioiain.machine.flows import predict



    data_path = None
    file = None
    chain = "A"

    if "--model" in sys.argv:
        data_path = sys.argv[sys.argv.index("--model") + 1]

    if "--file" in sys.argv:
        file = sys.argv[sys.argv.index("--file") + 1]

    if "--chain" in sys.argv:
        chain = sys.argv[sys.argv.index("--chain") + 1]

    print("MODEL PATH:", data_path)
    print("FILE PATH:", file)
    print("CHAIN:", chain)
    assert file is not None
    assert os.path.exists(file)

    if data_path is None:
        log("header", "Loading last model...")
        av_models = [os.path.join("./models",f) for f in os.listdir("./models") if f.endswith("data.json")]
        data_path = max(av_models, key = os.path.getctime)

    print("MODEL PATH:", data_path)
    print(f"Predicting contacts in file: {file}")

    SCRIPT = predict(file_path=file, model_data_path=data_path, chain_id=chain, use_temp=TEMP, timestamp=True)["script"]


if "-w" in sys.argv:
    log("title", "VISUALISATION")

    from src.bioiain.machine.flows import display_monomer_labels



    LABNAME = "dual_discrete_2"
    if "--label" in sys.argv:
        LABNAME = sys.argv[sys.argv.index("--label") + 1]





    if "--monomer" in sys.argv:
        target = sys.argv[sys.argv.index("--monomer") + 1]
        try:
            assert target in dataset
        except:
            [print(t) for t in dataset.embeddings.keys()]
            raise

        dataset.use_label(LABNAME)
        dataset.map()
        display_monomer_labels(monomer_id=target, dataset=dataset, label_name=LABNAME, script=SCRIPT, use_temp=TEMP)

    else:
        log("error", "Monomer id not provided (eg. --monomer 1M2Z_0_mon_A)")







log("end", "DONE")
log("title", "DONE")

end_pools()
sys.exit()