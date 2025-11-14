
from .crystal import Crystal, Path
from ..biopython import Model, Chain






class Oligomer(Model):
    child_class = Chain

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def _init(self, *args, **kwargs):
        pass


class OligomerBuilder(object):

    def build(self, crystal:Crystal, path:Path):
        print(crystal)
        print(path)
        self.path=path
