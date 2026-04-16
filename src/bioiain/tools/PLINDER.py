import os, sys
import pandas as pd

from ..utilities import *
from ..utilities.exceptions import *
from .. import TEMP_FOLDER, SUBDIR_NAME



class PLINDERSystem(object):
    def __init__(self, system_id=None):
        self.system_id = None
        self.system = None
        self.loaded = False
        if system_id is not None:
            self.from_id(system_id)

    def __repr__(self):
        return f"<PLINDERSystem: {self.system_id}>"

    def from_id(self, system_id):
        from plinder.core import PlinderSystem
        self.system_id = system_id
        self.system = PlinderSystem(system_id=system_id)
        self.loaded = True

    def _raise_if_unloaded(self):
        log("warning", "PLINDER system not loaded")
        raise PLINDERSystemNotLoaded()

    def annotations(self):
        self._raise_if_unloaded()
        return self.system.entry






class PLINDERDatabase(object):
    def __init__(self, tutorial=False):
        if tutorial:
            self.plinder_remote = "gs://plinder/2024-06/tutorial"
            self.plinder_mount = os.path.join(SUBDIR_NAME)
            self.plinder_cache = os.path.join(TEMP_FOLDER, "plinder/2024-06/tutorial")
            self.subset="tutorial"
        else:
            self.plinder_remote = "gs://plinder/2024-06/v2"
            self.plinder_mount = os.path.join(TEMP_FOLDER)
            self.plinder_cache = os.path.join(TEMP_FOLDER, "plinder/2024-06/v2")
            self.subset="v2"
        os.makedirs(self.plinder_cache, exist_ok=True)

        import plinder.core as pli
        self.pli = pli
        self.cfg = self.pli.get_config()

        self.cfg.data.plinder_mount = self.plinder_mount
        self.cfg.data.plinder_dir = self.plinder_cache
        self.cfg.data.plinder_remote = self.plinder_remote

        log("header", "Initialising PLINDER database")
        log(1, f"PLINDER remote data directory: {self.cfg.data.plinder_remote}")
        log(1, f"PLINDER mount directory: {self.cfg.data.plinder_mount}")
        log(1, f"PLINDER cache directory: {self.cfg.data.plinder_dir}")

    def __repr__(self):
        return f"<PLINDERDatabase ({self.subset}) at {self.plinder_mount}/plinder>"

    def query(self, columns=[], filters=[], **kwargs):
        return self.pli.scores.query_index(columns=columns, filters=filters, **kwargs)




















