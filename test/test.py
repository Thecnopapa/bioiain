import os, json, sys

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")

from src.bioiain.utilities.logging import *

tracemalloc_start()

from src.bioiain.aleph import *
from src.bioiain.base import *
from src.bioiain.machine import *
from src.bioiain.utilities.parallel import *

import torch, random
import numpy as np

torch.set_num_threads(avail_cpus)
log(1, f"Torch using {avail_cpus} threads")

seed = 6
random.seed(seed)
np.random.seed(seed)
torch.manual_seed(seed)
torch.cuda.manual_seed(seed)
torch.cuda.manual_seed_all(seed)

if "monomers" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "cath-monomeric",
                                  file_path="./data/cath-dataset-nonredundant-S20.monomeric.list",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/cath-monomeric"
    DATA_NAME = "monomers"
elif "receptors" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "receptors",
                                  file_path="./data/receptors.txt",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/receptors"
    DATA_NAME = "receptors"
elif "lbds" in sys.argv:
    if not "--no-download" in sys.argv:
        DATA_FOLDER = downloadPDB("./data", "lbds",
                                  file_path="./data/LBDs.txt",
                                  file_format="cif",
                                  overwrite=False)
    else:
        DATA_FOLDER = "./data/lbds"
    DATA_NAME = "lbds"
else:
    DATA_FOLDER = downloadPDB(list_name="aleph",
                              pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"],
                              data_dir="./data")
    DATA_NAME = "aleph"



for file in os.listdir(DATA_FOLDER):
    if "1M2Z" not in file.upper():
        continue

    entity = FragmentedStructure.from_file(os.path.join(DATA_FOLDER, file))

    matrix = entity.cvmatrix()
    entity.export()
    exit()
    entity.calculate_sasa()
    entity.export()
    embedding = CVEmbeddingV3C(entity=entity).embedding(force=True)
    entity.export()

    if embedding is None:
        log("warning", "No embedding for file:", file)
        continue

    [print(r, e) for r, e in zip(entity.residues(), embedding.tensor())]
        # entity.fragment().show_cvectors()






