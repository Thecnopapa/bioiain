import os, sys
from ..utilities import *
from .. import TEMP_FOLDER, SUBDIR_NAME


if "--tutorial" in sys.argv:
    plinder_remote = "gs://plinder/2024-06/tutorial"
    PLINDER_DIR = os.path.join(SUBDIR_NAME, "plinder/2024-06/tutorial")
else:
    plinder_remote = "gs://plinder/2024-06/v2"
    PLINDER_DIR = os.path.join(TEMP_FOLDER, "plinder/2024-06/v2")

import plinder.core
cfg = plinder.core.get_config()

cfg.data.plinder_dir = PLINDER_DIR
cfg.data.plinder_remote = plinder_remote

log("header", "Initialising PLINDER database")
log(1, f"PLINDER cache directory: {cfg.data.plinder_dir}")
log(1, f"PLINDER remote data directory: {cfg.data.plinder_remote}")
















