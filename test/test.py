import os, json, sys

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")

from src.bioiain.utilities.logging import *
from src.bioiain.base import *


if "monomers" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDBlist("./data", "cath-monomeric",
                                      file_path="./data/cath-dataset-nonredundant-S20.monomeric.list",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/cath-monomeric"
    DATA_NAME = "monomers"
elif "receptors" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDBlist("./data", "receptors",
                                      file_path="./data/receptors.txt",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/receptors"
    DATA_NAME = "receptors"
elif "lbds" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDBlist("./data", "lbds",
                                      file_path="./data/LBDs.txt",
                                      file_format="cif",
                                      overwrite=False)
    else:
        DATA_FOLDER = "./data/lbds"
    DATA_NAME = "lbds"
else:
    DATA_FOLDER = downloadPDBlist(list_name="aleph",
                                  pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"],
                                  data_dir="./data")
    DATA_NAME = "aleph"



from src.bioiain.aleph import *

entity = BIEntity.from_file(os.path.join(DATA_FOLDER, "7C2X.cif"), code="TEST", force=True)
print(entity)



exit()
from src.bioiain.machine.datasets import EmbeddingDataset

dataset = EmbeddingDataset(name=f"tokens_aleph_v4C").load()
print(dataset)
model = models.HopeLess(name="inference_test", in_shape=dataset.get(0).t.shape, inference=True)
model.load("/localdata/iain/bioiain/test/bioiain.d/models/HopeLess/HopeLess_monomers_v4C.temp.data.json")
print(model)

model.plot_latent_dimensions(dataset=dataset, r_threshold=100, only=[4,5], name="sasa_lig")


exit()



for file in os.listdir(DATA_FOLDER):
    if "1M2Z" not in file.upper():
        continue





    entity = FragmentedStructure.from_file(os.path.join(DATA_FOLDER, file))
    print(entity)
    entity.fragment()
    matrix = entity.cvmatrix()
    entity.export()
    entity.calculate_sasa()
    entity.export()
    from src.bioiain.machine import *
    embedding = CVEmbeddingV4C(entity=entity).embedding(force=True)
    entity.export()

    if embedding is None:
        log("warning", "No embedding for file:", file)
        continue
    print(embedding)
    for r, e in zip(entity.residues(), embedding.tensor()):
        print(r, "\t", " ".join([f"{ee.item():3.2f}" for ee in e]))
    #entity.show_cvectors()






