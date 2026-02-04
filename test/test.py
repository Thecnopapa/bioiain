import os, sys, json
sys.path.append('..')
from src.bioiain.biopython import downloadPDB
from src.bioiain import log


force = "force" in sys.argv

log("start", "SET UP")

#file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
#                          file_path="./pdb_list.txt", file_format="pdb",
#                          overwrite=False)
if "cath" in sys.argv:
    file_folder = downloadPDB("../../vib-ai/internship/data", "cath-nonredundant-S20",
                                           file_path="../../vib-ai/internship/data/cath-dataset-nonredundant-S20.list",
                                           file_format="cif",
                                           overwrite=False)
    pdb_list = "cath"
else:
    file_folder = downloadPDB("../../vib-ai/internship/data", "receptors",
                              file_path="../..//vib-ai/internship/data/receptors.txt", file_format="cif",
                              overwrite=False)
    pdb_list="rcps"

log("header", "File folder:", file_folder)


from src.bioiain.symmetries.elements import Monomer
from src.bioiain.symmetries.crystal import get_monomers
from src.bioiain.machine.datasets import EmbeddingDataset
from src.bioiain.machine.embeddings import SaProtEmbedding, MissingProgram, FoldseekError
from src.bioiain.symmetries.interactions import InteractionProfile
from src.bioiain.utilities.sequences import MSA


FORCE = "force" in sys.argv or "-f" in sys.argv

THRESHOLD=15
if "--threshold" in sys.argv:
    THRESHOLD = int(sys.argv[sys.argv.index("--threshold") + 1])

REBUILD = "--rebuild" in sys.argv
if FORCE:
    pass
ONLY = None
if "--only" in sys.argv:
    ONLY = sys.argv[sys.argv.index("--only") + 1].split(",")

data_name = f"saprot_interactions_{pdb_list}_T{THRESHOLD}"
log("header", f"Data name: {data_name}")

if ONLY is not None:
    data_name += "_ONLY_" + "_".join(ONLY)

dataset = EmbeddingDataset(name=data_name)
if not FORCE:
    dataset.load()




if "-l" in sys.argv or "-e" in sys.argv:

    if not "absolute_calcuated" in dataset.data or REBUILD:

        for n, file in enumerate(sorted(os.listdir(file_folder))):
            if not file.endswith(".cif"):
                continue

            if ONLY is not None:
                if file[:4] not in ONLY:
                    continue


            mon_data = get_monomers(file, file_folder, only_ids=True, force=FORCE, contact_threshold=15)
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
                    dataset.add_label_from_string(label, key=key)
                except MissingProgram as e:
                    raise e

                except FoldseekError as e:
                    log("warning", e)
                    continue
                except Exception as e:
                    log("Error", f"Exception occurred processing: {monomer_id}:\n", e)
                    raise e
                    continue

                dataset.save()


    dataset.data["absolute_calcuated"] = True
    dataset.save()


    datset_path = dataset.save()
    log("header", "DATASET:", dataset)
    msa = MSA(dataset.data["fasta_path"], dataset.data["name"])
    log("header", "MSA:", msa)

    if not "relative_calcuated" in dataset.data or REBUILD:

        for monomer_id in dataset.embeddings.keys():
            log("header", f"Calculating relative interactions for: {monomer_id}")

            if monomer_id in dataset:
                if "rel_label" in dataset.embeddings[monomer_id]:
                    if dataset.embeddings[monomer_id]["rel_label"] is not None and not (FORCE or REBUILD):
                        log(1, f"{monomer_id}: relative interactions already calculated")
                        continue
            monomer = Monomer.recover(data_path=os.path.join(dataset.data["export_folder"], monomer_id.split("_")[0], "monomers", monomer_id))
            if monomer is None:
                log("Warning", f"{monomer_id} has no monomer")
                continue
            log(1, "Generating relative labels...")
            ints = InteractionProfile(monomer, threshold=THRESHOLD, force=FORCE)
            rel_label = ints.generate_labels(relative=True, force=FORCE, dataset=dataset, msa=msa)
            dataset.add_label_from_list(rel_label, key=monomer_id, var_name="rel_label")
            print(dataset)


            dataset.save()

    dataset.data["relative_calcuated"] = True
    dataset.save()
    #if len(dataset) >= 1:
    #    log("start", "Dataset test")
    #    print(len(dataset))
    #    print(dataset[len(dataset)-1])





if "-t" in sys.argv:
    log("start", "Training")
    from src.bioiain.machine.models import *

    log("header", f"Dataset: {dataset}")
    dataset.use_label("rel_label")
    dataset.split()
    label_to_index = dataset.map(single_lab=True)
    log(2, "Label map:")
    print(json.dumps(label_to_index, indent=4))



    run_name = f"{data_name}"
    log(1, f"Run name: {run_name}")

    model = MLP_MK2(name=run_name, input_dim=1280, num_classes=len(label_to_index))
    model.add_map(dataset)

    epochs = 10
    if "--epochs" in sys.argv:
        epochs = int(sys.argv[sys.argv.index("--epochs") + 1])

    log("header", f"Epochs: {epochs}")


    model.add_epoch()
    for epoch in range(epochs):
        epoch = epoch +1
        log("header", "EPOCH:", epoch)
        dataset.train()
        try:
            for n, item in enumerate(dataset):
                #print("tensor:", item.t)
                #print("label:", item.l)

                out = model(item.t)
                #print("out:", out)

                loss = model.loss(out, item)
                model.step()


                print(f"{n:5d}/{len(dataset):5d}: LOSS={loss.item():.3f}", end = "\r")

            model.add_epoch()
            model.save()
            model.test(dataset)
        except KeyboardInterrupt:
            print("\nStopping model...")
            try:
                model.test(dataset)
            except ModelNotFound as e:
                print(e)
            break


if "-p" in sys.argv:
    from src.bioiain.symmetries import PredictedMonomerContacts
    from src.bioiain.visualisation import PymolScript

    prediction_folder = "./predictions"

    if "--file" in sys.argv:
        chains = None
        if "--chains" in sys.argv:
            chains = sys.argv[sys.argv.index("--chains") + 1]
        file = sys.argv[sys.argv.index("--file") + 1]
        print(f"Predicting contacts in file: {file}")
        assert os.path.exists(file)
        from src.bioiain.biopython import loadPDB
        structure = loadPDB(file, name=os.path.basename(file))
        print(structure)

        for chain in structure.get_chains():
            if chains is None or chain.id in chains:
                print(chain)
                monomer = Monomer.cast(chain)
                monomer.export(folder=prediction_folder)
                print(monomer)
                embedding = SaProtEmbedding(entity=monomer, folder=prediction_folder, force=FORCE)
                from src.bioiain.machine.models import *
                model = MLP_MK1(name="interactions", input_dim=480, num_classes=4)
                model.load("./models/MLP_MK2_saprot_interactions_rcps_T10.data.json")

                dataset = EmbeddingDataset(name=os.path.basename(file), folder=prediction_folder)
                dataset.add(embedding=embedding, key=monomer.get_name())
                label_to_index = model.data["label_to_index"]
                index_to_label = model.data["index_to_label"]
                print(dataset)
                full_pred=""
                with torch.no_grad():
                    for item in dataset:
                        out = model(item.t)
                        #print(out)
                        pred = out.max(dim=0)[1]
                        #print(pred)
                        p = index_to_label[str(pred.item())]
                        full_pred += p
                print(full_pred)
                interaction = PredictedMonomerContacts(monomer, full_pred[1:-1], label_to_index)
                pred_path = interaction.save_structure(prediction_folder)
                script = PymolScript(name=f"{monomer.get_name()}_prediction_pml_session", folder=prediction_folder)
                script.load(pred_path, monomer.get_name())
                script.spectrum(monomer.get_name())
                script.print(json.dumps(label_to_index, indent=4))
                script.write_script()
                #script.execute()
















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
exit()


