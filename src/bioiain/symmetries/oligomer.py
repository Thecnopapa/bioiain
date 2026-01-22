from copy import deepcopy

from . import entity_to_frac, coord_operation_entity, entity_to_orth, generate_displaced_copy
from .crystal import Crystal, Path
from ..biopython import Model, Chain
from ..utilities import log, add_front_0, find_com
from .operations import coord_add
from ..utilities import vector, rotation_matrix_from_vectors


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
        #if not(number == 3 and path["o_level"] == 4):
        #    return None
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
        self.chains.append([c for c in crystal.get_list() if c.id ==starting_chain][0].copy())
        #print(self.chains)
        if getattr(self.chains[0], "is_frac", False):
            self.chains[0] = entity_to_orth(self.chains[0], crystal.data["params"])


        # TODO: Work this out:
        z = list(zip(self.path["path"], self.path["coms"][2:]))
        for n, (step, com) in enumerate(z):
            print(n, step, com)
            chain_id = step["key"][-1]
            chain = [c for c in crystal.get_list() if c.id == chain_id][0].copy()
            chain.id = str(n)
            if not getattr(chain,"is_frac", False):
                entity_to_frac(chain, crystal.data["params"])

            old_vec = vector(crystal.data["symmetries"]["all_paths"][step["key"]]["vector"]["start"],
                             crystal.data["symmetries"]["all_paths"][step["key"]]["vector"]["end"])
            new_vec = vector(self.path["coms"][n+1], com)

            print("old_vec", old_vec)
            print("new_vec", new_vec)

            rot_mat = rotation_matrix_from_vectors(old_vec, new_vec)
            print(rot_mat)
            operation = {
                "rot": rot_mat,
                "tra": [0,0,0],
            }
            # coord_operation_entity(chain,
            #                        key=crystal.data["crystal"]["group_key"],
            #                        op_n=step["op_n"],
            #                        )
            coord_operation_entity(chain,operation=operation)



            # for sstep in self.path["path"][n:n]:
            #     coord_operation_entity(chain,
            #                            key=crystal.data["crystal"]["group_key"],
            #                            op_n=sstep["op_n"],
            #                            offset=sstep["position"],
            #                            )

            new_com = find_com(chain)
            delta = coord_add(com,new_com, True)
            generate_displaced_copy(chain, delta, copy=False)
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
        return oligomer




