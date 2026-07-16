
from ..utilities.logging import log

import torch, sys

DEVICE = "cpu"
DEVICE_N = None
USE_ALL_DEVICES = False
if "--cpu" not in sys.argv:
    if torch.cuda.is_available():
        DEVICE_LIST = list(range(torch.cuda.device_count()))
        DEVICE = "cuda"
        if "--cuda" in sys.argv:
            DEVICE_N =  [int(sys.argv[sys.argv.index("--cuda") + 1])]
            DEVICE = f"{DEVICE}:{DEVICE_N[0]}"
        elif "--cuda-all" in sys.argv and len(DEVICE_LIST) > 1:
            USE_ALL_DEVICES=True
            DEVICE_N = DEVICE_LIST
    elif torch.xpu.is_available():
        DEVICE = "xpu"

log(1, "DEVICE:", DEVICE)

FORCE = ("--force" in sys.argv) or ("-f" in sys.argv)
INFERENCE = ("--inference" in sys.argv) or ("-i" in sys.argv)
TRAIN = ("--train" in sys.argv) or ("-t" in sys.argv)
REBUILD = ("--rebuild" in sys.argv) or ("-r" in sys.argv)
EMBEDDINGS = ("--embeddings" in sys.argv) or ("-e" in sys.argv)
LABELS = ("--labels" in sys.argv) or ("-l" in sys.argv)


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
