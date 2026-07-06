
from ..utilities.logging import log

import torch, sys

DEVICE = "cpu"
if torch.cuda.is_available():
    DEVICE = "cuda"
    if "--cuda" in sys.argv:
        DEVICE_N =  int(sys.argv[sys.argv.index("--cuda") + 1])
        DEVICE = f"{DEVICE}:{DEVICE_N}"
elif torch.xpu.is_available():
    DEVICE = "xpu"

log(1, "DEVICE:", DEVICE)

FORCE = ("--force" in sys.argv) or ("-f" in sys.argv)
INFERENCE = "-i" in sys.argv
TRAIN = "-t" in sys.argv
REBUILD = ("--rebuild" in sys.argv) or ("-r" in sys.argv)


def set_seed(seed=6):
    import torch, numpy, random
    random.seed(seed)
    numpy.random.seed(seed)
    torch.manual_seed(seed)
    torch.cuda.manual_seed(seed)
    torch.cuda.manual_seed_all(seed)

def tensor_to_numpy(tensor):
    return tensor.detach().cpu().numpy()

from .embeddings import *
from .models import *
from .datasets import *
from .losses import *

__all_ = ["datasets", "embeddings", "losses", "models", "DEVICE"]
