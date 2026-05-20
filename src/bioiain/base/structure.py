from .entity import BIEntity
from .chain import BIChain





class BIStructure(BIEntity):
    child_class = BIChain
    extension = "structure"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def structure(self, *args, **kwargs):
        return self








