import os, json, sys


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


dataset = StructureDataset.from_list("./data/consensus.list")
print(dataset)
#dataset.load(entity_class=CompactStructure, check_existing=not FORCE)
#dataset.export()



DATASET_NAME = f"{dataset.name}_foldseek"
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
        print(t.shape)

        #print("DATA", data)
        #print("TENSOR", t.shape)
        #print("ENTITY:", entity)
        print("CHAIN:", chain)

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
            print("CHAIN 2:", chain)

            label_header="#"
            label_value=" "
            print(len(chain.residues()))

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

        print()


        total += 1


    print(f"Chains processed ok: {total-n_missmatches}, missmatches:{n_missmatches}")
    embeddings.save()

for e in embeddings:
    print(e)

print(embeddings)


exit()








model = CompactnessMLPmk1()
print(model)

epochs = 100
for epoch in range(epochs):
    print(f"EPOCH: {epoch}")


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



