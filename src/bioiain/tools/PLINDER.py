import os, sys
import pandas as pd

from ..utilities import *
from ..utilities.exceptions import *
from .. import TEMP_FOLDER, SUBDIR_NAME


if "--tutorial" in sys.argv:
    plinder_remote = "gs://plinder/2024-06/tutorial"
    PLINDER_MOUNT = os.path.join(SUBDIR_NAME)
    PLINDER_DIR = os.path.join(SUBDIR_NAME, "plinder/2024-06/tutorial")
else:
    plinder_remote = "gs://plinder/2024-06/v2"
    PLINDER_MOUNT = os.path.join(TEMP_FOLDER)
    PLINDER_DIR = os.path.join(TEMP_FOLDER, "plinder/2024-06/v2")
os.makedirs(PLINDER_DIR, exist_ok=True)

import plinder.core as pli
cfg = pli.get_config()

cfg.data.plinder_mount = PLINDER_MOUNT
cfg.data.plinder_dir = PLINDER_DIR
cfg.data.plinder_remote = plinder_remote

log("header", "Initialising PLINDER database")
log(1, f"PLINDER remote data directory: {cfg.data.plinder_remote}")
log(1, f"PLINDER mount directory: {cfg.data.plinder_mount}")
log(1, f"PLINDER cache directory: {cfg.data.plinder_dir}")



class PLINDERSystem(object):
    def __init__(self, system_id=None):
        self.system = None
        self.loaded = False
        if system_id is not None:
            self.from_id(system_id)

    @classmethod
    def from_id(self, system_id):
        from plinder.core import PlinderSystem
        self.system = PlinderSystem(system_id=system_id)
        self.loaded = True

    def _raise_if_unloaded(self):
        log("warning", "PLINDER system not loaded")
        raise PLINDERSystemNotLoaded()

    def annotations(self):
        self._raise_if_unloaded()
        return self.system.entry






class PLINDERDatabase(object):
    def __init__(self):
        pass

    def query(self, columns=[], filters=[], **kwargs):
        return pli.scores.query_index(columns=columns, filters=filters, **kwargs)




















