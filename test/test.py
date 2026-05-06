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

V0 = False
V1 = False
V2 = False
V3 = False
V4 = False

VA = False
VB = False
VC = False



if "--vB" in sys.argv:
    V0 = True
    VB = True
    DATA_NAME += "_vB"
    EMBEDDING_CLASS = CVEmbeddingVB

elif "--vC" in sys.argv:
    V0 = True
    VC = True
    DATA_NAME += "_vC"
    EMBEDDING_CLASS = CVEmbeddingVC

elif "--v1" in sys.argv:
    V1 = True
    VA = True
    DATA_NAME += "_v1"
    EMBEDDING_CLASS = CVEmbeddingV1
elif "--v2" in sys.argv:
    V2 = True
    VA = True
    DATA_NAME += "_v2"
    EMBEDDING_CLASS = CVEmbeddingV2

elif "--v1C" in sys.argv:
    V1 = True
    VC = True
    DATA_NAME += "_v1C"
    EMBEDDING_CLASS = CVEmbeddingV1C

elif"--v2C" in sys.argv:
    V2 = True
    VC = True
    DATA_NAME += "_v2C"
    EMBEDDING_CLASS = CVEmbeddingV2C

elif "--v3" in sys.argv or "--v3C" in sys.argv:
    V3 = True
    VC = True
    DATA_NAME += "_v3C"
    EMBEDDING_CLASS = CVEmbeddingV3C

else:
    DATA_NAME += "_v0"
    EMBEDDING_CLASS = CVEmbedding
    V0 = True
    VA = True

log(1, "DATA NAME:", DATA_NAME)



log("start", "Embeddings")
log("title", "Embeddings")

LR = 0.0001
if "--lr" in sys.argv:
    LR = float(sys.argv[sys.argv.index("--lr") + 1])
log(1, f"Learning rate: {LR}")

MODEL_NAME = "Hope"
if "--model" in sys.argv:
    MODEL_NAME = sys.argv[sys.argv.index("--model") + 1]
MODEL_CLASS = getattr(models, MODEL_NAME)
log(1, f"Model: {MODEL_CLASS}")

dataset = EmbeddingDataset(name=f"tokens_{DATA_NAME}")

if not ("--rebuild" in sys.argv or "--force" in sys.argv):
    dataset.load()
log(2, dataset)
total_files = len(os.listdir(DATA_FOLDER))
if len(dataset) == 0:
    log("ERROR", "Datest is empty:", dataset)
    exit()

#dataset.sequence_db()
#dataset.cluster(reassign=True, verbosity=3)
#dataset.save()

dataset.align(verbose=True, build_tree=True)
dataset.save()
log("end", "Embeddings")

model = None

log("start", "Tokenisation")
log("title", "Tokenisation")

if model is None:
    log("header", "Loading saved model...")
    model_data_path = sys.argv[sys.argv.index("--md") + 1]
    log(1, "Model path:", model_data_path)
    data = json.load(open(model_data_path))
    log(1, "Model class (data):", data.get("model"))
    model_class = getattr(models, data.get("model"))

    model = model_class(name="inference", in_shape=data.get("in_shape"), inference=True)
    model.load(model_data_path)
else:
    log("header", "Using loaded model...")

log(1, "Model:", model)
log(1, "Dataset:", dataset)

#tok_fasta = model._tokenise(dataset)



#matrix_path = model._build_blossum()
#model._align_tokens(dataset, tok_fasta, matrix="path", matrix_path=matrix_path, force=True)
model.plot_latent_space(show=True, mesh_points=50, dataset=dataset)









