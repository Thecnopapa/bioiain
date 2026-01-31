import os, sys
sys.path.append('..')



from src.bioiain.biopython import downloadPDB
from src.bioiain import log


force = "force" in sys.argv

log("start", "test.py")

#file_folder = downloadPDB("./data", "test_list", ["5JJM", "6nwl"],
#                          file_path="./pdb_list.txt", file_format="pdb",
#                          overwrite=False)
file_folder = downloadPDB("/home/iain/projects/vib-ai/internship/data", "receptors",
                          file_path="/home/iain/projects/vib-ai/internship/data/receptors.txt", file_format="cif",
                          overwrite=False)

log(1, "File folder:", file_folder)



from src.bioiain.symmetries.crystal import get_monomers
from src.bioiain.machine.datasets import EmbeddingDataset
from src.bioiain.machine.embeddings import SaProtEmbedding
from src.bioiain.symmetries.interactions import get_interaction_profile
dataset = EmbeddingDataset(name="saprot_with_interactions", tensor_iter_dim=1)


FORCE = "force" in sys.argv

for n, file in enumerate(sorted(os.listdir(file_folder))):
    if not file.endswith(".cif"):
        continue
    if "1M2Z" not in file:
        continue


    monomers = get_monomers(file, file_folder, force=FORCE)

    for monomer in monomers:
        print(">>>>", monomer)
        embedding = SaProtEmbedding(entity=monomer, force=FORCE)
        #print(embedding)
        key = dataset.add(embedding=embedding)
        label = get_interaction_profile(monomer, monomer.paths["export_folder"], threshold=10, force=FORCE)
        #print(key)
        #print(dataset.get(key, label=False))
        dataset.add_label_from_string(label, key=key)
        #print(dataset.get(key, label=True))


        print(dataset)

    log("start", "Dataset test")
    print(dataset[0])
    print(dataset[1])
    print(dataset[420])
    print(dataset)
    print("end")


    continue








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


