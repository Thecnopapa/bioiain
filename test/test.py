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

#file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
#                          file_path="./pdb_list.txt", file_format="pdb",
#                          overwrite=False)
if "no-download" in sys.argv:
    file_folder = "./data/cath-nonredundant-S20"
    pdb_list = "cath"
else:
    if "cath" in sys.argv:
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

DEVICE = "cpu"
if torch.cuda.is_available():
    DEVICE = "cuda"
elif torch.xpu.is_available():
    DEVICE = "xpu"

log(1, "DEVICE:", DEVICE)


FORCE = "force" in sys.argv or "-f" in sys.argv

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

data_name = f"saprot_interactions_{pdb_list}_T{THRESHOLD}"
if DUAL:
    data_name += "_DUAL"
    if DUAL_CLASSES:
        data_name += "_CLASSES"


if ONLY is not None:
    data_name += "_ONLY_" + "_".join(ONLY)

log("header", f"Data name: {data_name}")

dataset = EmbeddingDataset(name=data_name)
if not REBUILD:
    dataset.load()

from src.bioiain.machine.models import *

if "mk5" in sys.argv:
    model_class = DUAL_MLP_MK5
elif "mk6" in sys.argv:
    model_class = DUAL_MLP_MK6
else:
    if DUAL:
        model_class = DUAL_MLP_MK4
    else:
        model_class = MLP_MK3

log(1, f"Model class: {model_class}")

LR = 0.0001
if "--lr" in sys.argv:
    LR = float(sys.argv[sys.argv.index("--lr") + 1])

log(1, f"Learning rate: {LR}")

MIX = False
if "--mix" in sys.argv:
    MIX = True


#Comment

if "-l" in sys.argv or "-e" in sys.argv:

    if not "absolute_calcuated" in dataset.data or REBUILD:

        for n, file in enumerate(sorted(os.listdir(file_folder))):
            if not file.endswith(".cif"):
                continue

            if ONLY is not None:
                if file[:4] not in ONLY:
                    continue

            try:
                mon_data = get_monomers(file, file_folder, only_ids=True, force=FORCE, contact_threshold=15)
            except (FileNotFoundError, StructureRecoverException):
                print("retrying")
                mon_data = get_monomers(file, file_folder, only_ids=True, force=True, contact_threshold=15)
                monomer = Monomer.recover(data_path=os.path.join(dataset.data["export_folder"], monomer_id.split("_")[0], "monomers", monomer_id))

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
                    continue

                dataset.save()


    dataset.data["absolute_calcuated"] = True
    dataset.save()
    if DUAL_CLASSES:
        dataset.use_label("abs_label")
        dataset.map()


    datset_path = dataset.save()
    log("header", "DATASET:", dataset)
    msa = MSA(dataset.data["fasta_path"], dataset.data["name"])
    log("header", "MSA:", msa)

    if (not "relative_calcuated" in dataset.data) or REBUILD or "-rel" in sys.argv:

        for monomer_id in list(dataset.embeddings.keys()):
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

                rel_label = ints.generate_labels(relative=True, force=FORCE, dataset=dataset, msa=msa, dual=DUAL, in_lab_var="abs_label")
                print("REL LAB:", len(rel_label), f"DUAL={DUAL}")

                if DUAL:
                    if DUAL_CLASSES:
                        dataset.add_label_from_list(rel_label, key=monomer_id, var_name="dual_class_label")
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


    if DUAL_CLASSES:
        dataset.use_label("dual_class_label")
    elif DUAL:
        dataset.use_label("dual_label")
    else:
        dataset.use_label("rel_label")

    dataset.map()
    dataset.data["relative_calcuated"] = True
    dataset.save()
    #if len(dataset) >= 1:
    #    log("start", "Dataset test")
    #    print(len(dataset))
    #    print(dataset[len(dataset)-1])





if "-t" in sys.argv:
    log("title", "TRAINING")
    log("start", "Training")

    log("header", f"Dataset: {dataset}")
    if DUAL_CLASSES:
        dataset.use_label("dual_class_label")
    elif DUAL:
        dataset.use_label("dual_label")
    else:
        dataset.use_label("rel_label")
    dataset.split()
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

    
    run_name = f"{data_name}_LR{str(LR).split(".")[-1]}_E{epochs}_MIX:{int(MIX)}"
    model = model_class(name=run_name, in_shape=(1280,), num_classes=len(label_to_index), lr=LR, weights=label_count  ).to(DEVICE)
    model.add_map(dataset)
    run_name = model.data["name"]


    log(1, f"Run name: {run_name}")



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

        if not is_cluster:
            epoch_start = datetime.datetime.now()
            eta = datetime.timedelta(0)
        try:
            for n, item in enumerate(dataset):
                if not is_cluster:
                    print(f"\033]0;TRAIN {epoch}/{epochs} {(n/len(dataset))*100:3.0f}%\a", end="\r")
                #print("tensor:", item.t)
                #print("label:", item.l)
                print(DEVICE)
                item.to(DEVICE)
                out = model(item.t)
                print("out:", out)
                print(type(item.l), item.l)

                loss = model.loss(out, item)
                #print(out.item(), item.l, loss)
                #exit()
                model.step()

                if not is_cluster:
                    if n % (data_len // 100) == 0:
                        now = datetime.datetime.now() 
                        delta = now - epoch_start
                        eta = datetime.timedelta(seconds=(delta.total_seconds() / (n+1)) * (data_len - n))
                    print(f"{n:6d}/{len(dataset):6d}: LOSS={model.running_loss['default']/model.running_loss['total']:6.4f} ({loss:6.4f}) <-- o:{torch.max(out, dim=0)[1].item()} t:{item.li} ETA: {eta}", end = "\r")

            model.save()
            model.test(dataset)
            model.send_run(host="iainvisa.com", key=os.environ.get("IAINVISA_FILE_KEY"), epoch=epoch)
            model.add_epoch()

        except KeyboardInterrupt:
            log("warning", "!! Stopping model...")
            try:
                #print(os.environ)
                #print("KEY:", os.getenv("IAINVISA_FILE_KEY"))
                #model.send_run(host="iainvisa.com", key=os.environ.get("IAINVISA_FILE_KEY"))
                #print(model.send_run(host="127.0.0.1:5000", key=os.getenv("IAINVISA_FILE_KEY"), protocol="http").__dict__)
                print(repr(model))
                model.test(dataset, re_load=False)
            except ModelNotFound as e:
                print(e)
            break


if "-p" in sys.argv:
    from src.bioiain.symmetries import PredictedMonomerContacts
    from src.bioiain.visualisation import PymolScript
    log("title", "PREDICTING...")
    prediction_folder = "./predictions"

    model_path = "./models/DUAL_MLP_MK3_saprot_interactions_rcps_T10_DUAL_CLASSES.temp.data.json"
    if "--model" in sys.argv:
        model_path = sys.argv[sys.argv.index("--model") + 1]

    if "--file" in sys.argv:
        chains = None
        if "--chains" in sys.argv:
            chains = sys.argv[sys.argv.index("--chains") + 1]
        file = sys.argv[sys.argv.index("--file") + 1]
        print(f"Predicting contacts in file: {file}")
        assert os.path.exists(file)
        from src.bioiain.biopython import loadPDB
        structure = loadPDB(file, name=os.path.basename(file).split(".")[0])
        print(structure)

        for chain in structure.get_chains():
            if chains is None or chain.id in chains:
                print(chain)
                monomer = Monomer.cast(chain)
                monomer.export(folder=prediction_folder)
                print(monomer)
                embedding = SaProtEmbedding(entity=monomer, folder=prediction_folder, force=True)
                from src.bioiain.machine.models import *
                data = json.load(open(model_path))
                weights = data.get("weights", [])
                model = model_class(name="interactions", in_shape=data["in_shape"], num_classes=data["num_classes"], weights=weights)
                del data
                model.load(model_path)
                print(model)
                print(json.dumps(model.data, indent=4))

                dataset = EmbeddingDataset(name=os.path.basename(file), folder=prediction_folder)
                dataset.add(embedding=embedding, key=monomer.get_name())
                label_to_index = model.data["label_to_index"]
                index_to_label = model.data["index_to_label"]
                print(dataset)

                full_pred=[]
                with torch.no_grad():
                    for item in dataset:
                        out = model(item.t)
                        #print(out)
                        #print(out.shape)
                        if out.shape[0]>2:
                            pred = out.max(dim=0)[1]
                        else:
                            pred = out
                        #print(pred)
                        if len(label_to_index) > 2:
                            p = index_to_label[str(pred.item())]
                            full_pred.append(p)
                        elif len(label_to_index) == 2:
                            cp = round(pred[0].item()*100)
                            op = int(pred[1].item() > 0.5)
                            full_pred.append(pred[1].item())
                        else:
                            full_pred.append(pred.item())
                #print(full_pred)

                interaction = PredictedMonomerContacts(monomer, full_pred, label_to_index)
                pred_path = interaction.save_structure(prediction_folder)
                script = PymolScript(name=f"{monomer.get_name()}_{chain.id}_prediction_pml_session", folder=prediction_folder)
                script.load(pred_path, monomer.get_name())
                script.spectrum(monomer.get_name())
                if label_to_index is not None:
                    script.print(json.dumps(label_to_index, indent=4))
                session_path = script.write_script()
                print("Session saved at:")
                print("pymol", session_path)
                #script.execute()

if "-w" in sys.argv:
    from src.bioiain.symmetries import PredictedMonomerContacts
    from src.bioiain.visualisation import PymolScript

    LABNAME = "dual_class_label"
    if "--label" in sys.argv:
        LABNAME = sys.argv[sys.argv.index("--label") + 1]

    dataset.use_label(LABNAME)
    dataset.map()


    if "--monomer" in sys.argv:
        target = sys.argv[sys.argv.index("--monomer") + 1]
        try:
            assert target in dataset
        except:
            [print(t) for t in dataset.embeddings.keys()]
            raise
        log("title", "VISUALISATION")
        log("start", "VISUALISATION")
        log("header", f"Displaying monomer: {target}")

        embedding = dataset.embeddings[target]
        print(json.dumps(embedding, indent=4))

        label_path = embedding[LABNAME]
        with open(label_path, "r") as f:
            label = f.read().split(",")

        monomer = Monomer.recover(data_path=os.path.join(dataset.data["export_folder"], target.split("_")[0], "monomers", target))


        log(1, f"label: {label}")
        print(dataset.data["label_to_index"])
        if len(label) == 1:
            label = label[0]
        interaction = PredictedMonomerContacts(monomer, label, label_to_index=dataset.data["label_to_index"])
        view_path = interaction.save_structure("/tmp/bioiain/visualisations", extra_name="_visualised_monomer_contacts")
        log(1, interaction)

        script = PymolScript(name=f"{monomer.get_name()}_visualisation_pml_session", folder="/tmp/bioiain/visualisations")
        script.load(view_path, monomer.get_name())
        script.spectrum(monomer.get_name())


        session_path = script.write_script()
        print("Session saved at:")
        print("pymol", session_path)














#     script.line(f"int_{n}", sele1=f"monomer and c. {a1['chain']} and i. {a1['resn']} and n. CA",
#                 sele2=f"interacting_{n} and c. {a2['chain']} and i. {a2['resn']} and n. CA")
# script.disable(name)

# script = pymol.PymolScript(name=monomer, pymol_path="$CONDA_PREFIX/bin/pymol")
# script.load(crystal.paths["original"], "original", to="pdb")
# script.cell()
# script.symmetries()
# script.group()
# script.disable("sym")
# script.disable("original")


# embeddings.append(interactions_per_monomer(monomer, crystal.paths["monomer_folder"], script=script))

# script.write_script()
# script.execute()




# Oligo  (to be reworked)


    # crystal.get_oligomers(
    #     oligomer_levels=[2],
    # )


    # from src.bioiain.symmetries import Oligomer
    #
    # print(crystal.paths["oligo_folder"])
    # for file in sorted(os.listdir(crystal.paths["oligo_folder"]))   :
    #     print(file)
    #     if file.endswith(".data.json"):
    #         oligo = Oligomer.recover(id="recovered",data_path=os.path.join(crystal.paths["oligo_folder"], file))
    #         print(oligo)
    #
    #
    # from src.bioiain.visualisation import pymol
    #
    # script = pymol.PymolScript(folder=".", name="test")
    # script.load(crystal.paths["original"], "original", to="pdb")
    # script.cell()
    # script.symmetries()
    # script.group()
    # script.write_script()






log("end", "DONE")
log("title", "DONE")

end_pools()
sys.exit()