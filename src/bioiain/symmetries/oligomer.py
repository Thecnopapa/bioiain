from copy import deepcopy

from . import entity_to_frac, coord_operation_entity, entity_to_orth
from .crystal import Crystal, Path
from ..biopython import Model, Chain
from ..utilities import log, add_front_0


class Oligomer(Model):
    child_class = Chain

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def _init(self, *args, **kwargs):
        self.paths["export_folder"] = self.paths["export_folder"].replace("crystal", "oligomers")
        self.data["info"]["name"] = self.data["info"]["name"].replace("cryst", "oligo")



class OligomerBuilder(object):

    def build(self, crystal:Crystal, path:Path, number:int):
        log(3, "Building Oligomer number {}".format(number))
        log(4, crystal)
        log(4, path)
        self.path=path
        self.model=Model(number)
        self.model.base_init()
        self.model.data["info"] = deepcopy(crystal.data["info"])
        self.model.paths = deepcopy(crystal.paths)
        log(4, self.model)
        self.chains = []
        starting_chain = self.path["path"][0]["key"][0]
        log(5, "starting chain:", starting_chain)
        self.chains.append([c.copy() for c in crystal.get_list() if c.id ==starting_chain][0])
        #print(self.chains)
        if getattr(self.chains[0], "is_frac", False):
            self.chains[0] = entity_to_orth(self.chains[0], crystal.data["params"])


        # TODO: Work this out:
        for n, step in enumerate(self.path["path"]):
            chain_id = step["key"][-1]
            chain = [c.copy() for c in crystal.get_list() if c.id == chain_id][0]
            chain.id = str(n)
            if not getattr(chain,"is_frac", False):
                entity_to_frac(chain, crystal.data["params"])
            for sstep in self.path["path"][:n]:
                coord_operation_entity(chain,
                                       key=crystal.data["crystal"]["group_key"],
                                       op_n=sstep["op_n"]
                                       )

            entity_to_orth(chain, crystal.data["params"])
            self.chains.append(chain)

        for chain in self.chains:
            self.model.add(chain)

        oligomer = Oligomer.cast(self.model)
        oligomer.data["info"]["name"] = oligomer.data["info"]["name"] + "_{}_{}".format(self.path["o_level"], add_front_0(number, digits=3))
        oligomer.data["path"] = path
        path_keys = [p["key"] for p in path["path"]]
        oligomer.data["paths"] = {k:v for k,v in crystal.data["symmetries"]["all_paths"].items() if k in path_keys}
        oligomer.data["params"] = crystal.data["params"]
        oligomer.data["crystal"] = crystal.data["crystal"]


        oligomer.pass_down()
        oligomer.export()
        log(4, oligomer)
        return oligomer.paths["export_folder"]




