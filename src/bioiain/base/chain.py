import os, json

from ..utilities.exceptions import *
from .entity import BIEntity
from .residue import BIResidue

class BIChain(BIEntity):
    child_class = BIResidue
    extension = "chain"
    level = "chain"

    def __init__(self):
        super().__init__()
        self.data["export_subdir"] = "chains"