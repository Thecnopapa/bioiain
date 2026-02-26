__all_ = ["datasets", "embeddings", "models"]

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