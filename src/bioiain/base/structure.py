from .entity import Entity
from .chain import Chain





class Structure(Entity):
    child_class = Chain
    extension = "structure"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def structure(self, *args, **kwargs):
        return self








