import os, sys, json

from ..utilities.exceptions import *
from .entity import BIEntity
from .chain import BIChain





class BIStructure(BIEntity):
    child_class = BIChain
    extension = "structure"

    def __init__(self):
        super().__init__()




