import os, sys, json
from textwrap import indent

sys.path.append('..')



from src.bioiain.biopython import downloadPDB
from src.bioiain import log


force = "force" in sys.argv

log("start", "test.py")

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

log(1, "File folder:", file_folder)


from src.bioiain.symmetries.elements import Monomer
from src.bioiain.symmetries.crystal import get_monomers
from src.bioiain.machine.datasets import EmbeddingDataset
from src.bioiain.machine.embeddings import SaProtEmbedding, MissingProgram
from src.bioiain.symmetries.interactions import get_interaction_profile


FORCE = "force" in sys.argv or "-f" in sys.argv
if FORCE:
    pass

run_name = f"saprot_interactions_{pdb_list}"

dataset = EmbeddingDataset(name=run_name)
if not FORCE:
    dataset.load()


if "-l" in sys.argv or "-e" in sys.argv:

    for n, file in enumerate(sorted(os.listdir(file_folder))):
        if not file.endswith(".cif"):
            continue
        #if "1M2Z" not in file:
        #    continue


        mon_data = get_monomers(file, file_folder, only_ids=True, force=FORCE)

        print(mon_data)
        if mon_data is None:
            log("Warning", f"{file} has no monomers")
            continue

        monomers, monomer_folder = mon_data

        for monomer_id in monomers:
            try:
                print(">>>>", monomer_id)
                if monomer_id in dataset and not force:
                    print(f"{monomer_id} already in dataset")
                    continue

                monomer = Monomer.recover(data_path=os.path.join(monomer_folder, monomer_id))
                if monomer is None:
                    log("Warning", f"{monomer_id} has no monomer")
                    continue
                embedding = SaProtEmbedding(entity=monomer, force=FORCE)
                key = dataset.add(embedding=embedding, key=monomer.get_name())
                label = get_interaction_profile(monomer, monomer.paths["export_folder"], threshold=10, force=FORCE)
                dataset.add_label_from_string(label, key=key)
                print(dataset)
            except MissingProgram as e:
                raise e
            except Exception as e:
                log("Error", f"Exception occurred processing: {monomer_id}:\n", e)
                #raise e
                continue

            dataset.save()

        if len(dataset) >= 1:
            log("start", "Dataset test")
            print(len(dataset))
            print(dataset[len(dataset)-1])

        continue

    datset_path = dataset.save()
    print(dataset, f"saved at: {datset_path}")


elif "-t" in sys.argv:
    print(dataset)
    dataset.split()

    label_to_index = dataset.map()
    print(json.dumps(label_to_index, indent=4))


    # dataset.train()
    # import random, math
    # for n in range(10):
    #     key = math.floor(random.random() * len(dataset))
    #     print("KEY:", key)
    #     print(dataset[key])
    #     print()


    from torch.utils.data import DataLoader
    import torch.nn as nn
    import torch.optim as optim

    dataloader = dataset

    from src.bioiain.machine.models import *
    model = MLP_MK1(name=run_name, input_dim=480, num_classes=len(label_to_index))
    model.add_map(dataset)


    criterion = nn.CrossEntropyLoss()
    optimizer = optim.Adam(model.model.parameters())

    epochs = 10







    log("start", "Training")
    for epoch in range(epochs):
        log("header", "EPOCH:", epoch)
        dataset.train()
        try:
            for n, item in enumerate(dataset):
                #print("tensor:", item.t)
                #print("label:", item.l)

                truth = [0] * len(label_to_index)
                truth[label_to_index[item.l]] = 1

                out = model(item.t)
                #print("out:", out)

                optimizer.zero_grad()
                loss = criterion(out, torch.Tensor(truth))

                loss.backward()
                optimizer.step()
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


elif "-p" in sys.argv:

    if "--file" in sys.argv:
        chains = None
        if "--chains" in sys.argv:
            chains = sys.argv[sys.argv.index("--chains") + 1]
        file = sys.argv[sys.argv.index("--file") + 1]
        print(f"Predicting contacts in file: {file}")
        assert os.path.exists(file)
        from src.bioiain.biopython import loadPDB
        structure = loadPDB(file, name="prediction",)
        print(structure)

        for chain in structure.get_chains():
            if chains is None or chain.id in chains:
                print(chain)
                monomer = Monomer.cast(chain)
                monomer.export()
                print(monomer)
                embedding = SaProtEmbedding(entity=monomer, force=FORCE)
                from src.bioiain.machine.models import *
                model = MLP_MK1(name="interactions", input_dim=480, num_classes=4)
                model.load("./models/MLP_MK1_interactions.data.json")

                dataset = EmbeddingDataset(name="prediction")
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


