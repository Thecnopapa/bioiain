
from .crystal import Crystal, Path
from ..biopython import Model, Chain
from ..utilities import log






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
        log(3, "Building Oligomer")
        log(4, crystal)
        log(4, path)
        self.path=path
        self.model=Model(0)
        log(4, self.model)
        self.chains = []
        starting_chain = self.path["path"][0]["key"][0]
        log(5, "starting chain:", starting_chain)
        self.chains.append([c for c in crystal.get_list() if c.id ==starting_chain][0])
        print(self.chains)
