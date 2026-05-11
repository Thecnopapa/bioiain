
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


def tensor_to_numpy(tensor):
    return tensor.detach().cpu().numpy()

from .embeddings import *
from .base_model import *
from .datasets import *
from .flows import *
from .losses import *
from .models import *

__all_ = ["datasets", "embeddings", "models", "flows", "losses", "base_model"]
