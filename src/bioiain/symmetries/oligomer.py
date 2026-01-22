from copy import deepcopy

from . import entity_to_frac, coord_operation_entity, entity_to_orth, generate_displaced_copy
from .crystal import Crystal
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




def _calculate_oligomerisation_paths(self, show=False) -> Self:
        if self.monomers is None:
            self._identyfy_main_elements()
            self._regenerate_crystal()

        log(1, "Calculating oligomerisation paths ({})".format(self.data["info"]["name"]))

        fig, ax = self.plot(show=False)

        for n, monomer in enumerate(self.monomers):
            log(2, "Monomer: {}".format(monomer.data["info"]["name"]))
            for oc in monomer.data["contacts"]["all"]:
                log(3, oc["name"])
                com1 = oc["monomer1"]["CoM-frac"]
                com2 = oc["monomer2"]["CoM-frac"]
                if monomer.data["contacts"]["paths"] == None:
                    monomer.data["contacts"]["paths"] = {}

                if com2 is None:
                    if oc["monomer2"]["operation"] != None:
                        target_mon = [m for m in self.monomers if m.id == oc["monomer2"]["id"]][0]
                        for position in oc["positions"]:
                            c = deepcopy(oc)
                            c["name"] = c["name"] + ", pos_{}".format("".join([str(p) for p in position]))
                            c["position"] = position
                            log(4, "Monomer 2 at:", position)
                            displaced_monomer = generate_displaced_copy(
                                target_mon.copy(),
                                distance=position,
                                key=self.data["crystal"]["group_key"],
                                op_n=c["monomer2"]["operation"])
                            com2 = find_com(displaced_monomer)
                            c["vector"] = {"start": com1, "end": com2}
                            ax.add_artist(Arrow3D(*zip(com1, com2), color=pymol_colours[n]))
                            monomer.data["contacts"]["paths"]["{}{}{}".format(
                                monomer.id,
                                len(monomer.data["contacts"]["paths"].keys()),
                                displaced_monomer.id)] = c
                    else:
                        c = deepcopy(oc)
                        c.name = c["name"] + "pos_{}".format("AU")
                        c["position"] = "AU"
                        com2 = find_com([m for m in self.monomers if m.id == c["monomer2"]["id"]][0])
                        log(4, "Monomer 2 at:", c["positions"], "(Asymmetric unit)")
                        ax.add_artist(Arrow3D(*zip(com1, com2),color=pymol_colours[n]))
                        monomer.data["contacts"]["paths"]["{}{}".format(
                            monomer.id,
                            len(monomer.data["contacts"]["paths"].keys()))] = c
        self.data["symmetries"]["all_paths"] = {}
        for m in self.monomers:
            if m.data["contacts"]["paths"] is not None:
                #print(m.data.keys())
                #print(m.data["contacts"])
                self.data["symmetries"]["all_paths"].update(m.data["contacts"]["paths"])
        self.export()
        if show:
            fig.show()
            input("Press Enter to continue...")


        return self


def _find_oligomers(self, plot=False) -> Self:
        if not self.data["symmetries"].get("all_paths", False):
            self._calculate_oligomerisation_paths()

        log(1, "Pathing oligomers ({})".format(self.data["info"]["name"]))
        log(2, "Oligomer levels: {}".format(self.data["crystal"]["oligomer_levels"]))

        log(2, "Available paths:")
        for k, contact in self.data["symmetries"]["all_paths"].items():
            log(3,k, ":", contact["name"])
        total_paths = []

        p = Path(self.data["info"]["name"])
        total_paths.extend(self._extend_path(p,
                                             self.data["symmetries"]["all_paths"],
                                             self.data["crystal"]["oligomer_levels"]))
        log(2, "Total paths:", len(total_paths))

        unique_paths = []

        for path in total_paths:
            steps = [k["key"] for k in path["path"]]
            path["steps"] = steps
            exists = False
            if steps in unique_paths or steps.copy().reverse() in unique_paths:
                exists = True
                continue
            for i in range(1,len(steps)-1):
                if steps[i:]+steps[:i] in unique_paths:
                    exists = True
            if not exists:
                unique_paths.append(path)


        log(2, "Unique paths:", len(unique_paths))
        [log(3, path["steps"]) for path in unique_paths]
        log(2, "Unique paths:", len(unique_paths))

        all_paths = self.data["symmetries"]["all_paths"]
        print(all_paths.keys())
        print(self.monomers)
        mons = {
            m.id: m for m in self.monomers
        }
        o_coms = {
            m.id: find_com(m) for m in self.monomers
        }

        final_paths = {k:[] for k in self.data["crystal"]["oligomer_levels"]}

        for n, path in enumerate(unique_paths):
            omit = False
            relevant = False

            log(3, "Path: {}, level: {}".format(n, path["o_level"]) )
            log(4, "Steps: {}".format(path["steps"]))
            starting_monomer = path["steps"][0][0]
            log(4, "Starting monomer: {}".format(starting_monomer))

            if plot:
                fig, ax = self.plot(show=False)

            point_list = []
            operation_list = []
            current_pos = [0, 0, 0]

            coms = deepcopy(o_coms)
            #print("COMS:", coms)
            point_list.append(coms[starting_monomer])

            for s, step in enumerate(path["path"]):
                step_info = all_paths[step["key"]]
                #print("###")
                #print(step_info)
                #print("###")
                #print(step)
                #print("###")
                op_n = step_info["monomer2"]["operation"]
                key = self.data["crystal"]["group_key"]
                params = self.data["params"]
                pos = step_info["position"]

                reverse = step["reverse"]
                # Fix reverse / Revisit revesed / Skipped for now
                if reverse:
                    log("warning", "Reverse operations not yet supported, will revisit if needed ")
                    omit = True
                    break

                #print("pos:", pos)
                #print("op_n:", op_n)
                #print("op_list:", operation_list)

                id1 = step_info["monomer1"]["id"]
                id2 = step_info["monomer2"]["id"]



                origin = [0,0,0]

                if len(operation_list) == 0:
                    if plot:
                        ax.scatter(*origin, color="purple")
                        ax.text(*origin, "o", color="purple")
                    com1 = step_info["vector"]["start"]
                    com2 = step_info["vector"]["end"]
                    position = com2
                    point_list.append(com1)
                    point_list.append(com2)

                else:
                    com1 = operation_list[-1]["position"]

                    vec = deepcopy(step_info["vector"])
                    for o, op in enumerate(operation_list):
                        for k, v in vec.items():
                            vec[k] = coord_operation(
                                v,
                                key=key,
                                op_n=op["op_n"],
                                #distance=op["pos"],
                            )

                    if reverse:
                        for k, v in vec.items():
                            vec[k] = coord_operation(
                                v,
                                key=key,
                                op_n=op["op_n"],
                                #distance=op["pos"],
                                reverse=True
                            )

                    if not reverse:
                        c2 = vector(vec["start"], vec["end"])
                    else:
                        c2 = vector(vec["end"], vec["start"])
                    #print("c2:", c2)



                    com2 = coord_add(com1, c2)
                    position = com2
                    ca = "".join([str(round(c, 3)) for c in com2])
                    cl = ["".join([str(round(c, 3)) for c in cp]) for cp in point_list]
                    print(ca)
                    print(cl)
                    if ca in cl:
                        omit = True
                        log("warning", "Path does walkback".format(path["steps"]))
                        break

                    point_list.append(com2)

                if omit:
                    log("warning", "Path: {} discarded".format(path["steps"]))
                    break


                #print("Com1:", com1)
                #print("Com2:", com2)

                if reverse:
                    #print(id1, "<--", id2)
                    pass
                else:
                    #print(id1, "-->", id2)
                    pass


                operation_list.append({
                    "op_n": op_n,
                    "pos": pos,
                    "key": key,
                    "vector": step_info["vector"],
                    "position": position,
                    "reverse": reverse,

                })

                if plot:
                    if reverse:
                        ax.add_artist(Arrow3D(*zip(com1, com2),color="blue", alpha=0.5))
                        ax.add_artist(Arrow3D(*zip(*step_info["vector"].values()),color="green", alpha=0.5))
                    else:
                        ax.add_artist(Arrow3D(*zip(com1, com2),color="red", alpha=0.5))


                # Check for early circle closures -> Remove if so

            # Then classify linear/circular based on final point (whether is on direct contact with ASU)
            # -> final point in vector ends list
            relevant = True

            if (not omit) and relevant:
                log(4, "Saving path...")
                path["coms"] = point_list
                if plot:
                    fig_path = os.path.join(self.paths["export_folder"], "figs")
                    os.makedirs(fig_path, exist_ok=True)
                    path["fig_path"] = fig.savefig(os.path.join(fig_path,"{}_path_{}_{}.png".format(
                        self.data["info"]["name"],
                        path["o_level"], n)))
                    log(4, "Figure saved at: {}".format(path["fig_path"]))
                final_paths[path["o_level"]].append(path)


        self.data["symmetries"]["relevant_paths"] = final_paths
        log(2, "Relevant paths:")
        for ol in final_paths.keys():
            log(3, "Oligomer level:", ol)
            log(4, "Total:", len(final_paths[ol]))
            [log(5, path["steps"]) for path in final_paths[ol]]
            log(4, "Total:", len(final_paths[ol]))
        return self


def _extend_path(path, options, o_levels, depth=3):
        branches = []
        o_levels = deepcopy(o_levels)
        if path.length in o_levels:
            path = deepcopy(path)
            path.close()
            branches.append(deepcopy(path.data))
            o_levels.remove(path.length)
            if len(o_levels) == 0:
                log(depth, "DONE", path.length, o_levels, len(o_levels))
                return branches
        log(depth, path.length, o_levels, len(o_levels))


        last_id = None
        if path.length == 1:
            available_options = options
        else:
            last_path = path.data["path"][-1]
            last_key = last_path["key"]
            if last_path["reverse"]:
                last_id = last_key[0]
            else:
                last_id = last_key[-1]

            available_options = {k:o for k, o in options.items()
                                 if (k.startswith(last_id)
                                     #or k.endswith(last_id) # Reverse disabled
                                     ) and k != last_key}
            #log(depth,"Available", available_options.keys(), last_path, last_mon)



        reverse = False
        for k, v in available_options.items():
            log(depth, "Branching to {}".format(k))
            branch_path = deepcopy(path)
            if last_id is not None:
                reverse = not k.startswith(last_id)
            branch_path.add(k,v, reverse=reverse)
            branches.extend(Crystal._extend_path(branch_path, options, o_levels, depth=depth+1))


        return branches

def _build_oligomers(self):
        log(1, "Building oligomers ({})".format(self.data["info"]["name"]))
        oligo_folder = None
        from .oligomer import OligomerBuilder
        if not self.data["symmetries"].get("relevant_paths", False):
            self._find_oligomers()
        for ol in self.data["symmetries"]["relevant_paths"].keys():
            log(2, "O level:", ol)
            for n, path in enumerate(self.data["symmetries"]["relevant_paths"][ol]):
                builder = OligomerBuilder()
                oligo = builder.build(self, path, number=n)
                self.paths["oligo_folder"] = oligo.paths["export_folder"]
        return self




class Path(object):
    def __init__(self, name):
        self.data = {}
        self.data["name"] = name
        self.data["path"] = []
        self.contacts = []
        self.data["o_level"] = None
        self.data["complete"] = False
        self.length = 1
        self.complete = False


    def close(self):
        self.complete = True
        self.data["complete"] = True
        self.data["o_level"] = len(self.data["path"])+1

    def add(self, path, contact, reverse=False ):
        self.length += 1
        self.data["path"].append({
            "key": path,
            "reverse":reverse,
            "op_n": contact["monomer2"]["operation"],
            "pos": contact["position"],
        })
        self.contacts.append(contact)

    def __repr__(self):
        return "<Path in {}: {}, level: {}, closed: {}>".format(
            self.data["name"],
            "-".join([k["key"] for k in self.data["path"]]),
            self.data["o_level"],
            self.data["complete"])
