
import os


from typing_extensions import Self
from copy import deepcopy

from .operations import *
from ..utilities.logging import log
from ..utilities.maths import *
from ..biopython import Model

from ..visualisation import fig3D, pymol_colours, Arrow3D



class Crystal(Model):
    def _init(self, *args, **kwargs) -> Self:
        """
        Initialize the crystal. Executed on casting.
        :param args:
        :param kwargs:
        :return: Self.
        """
        super()._init(*args, **kwargs)
        self.data["crystal"]["oligomer_levels"] = None
        self.data["crystal"]["min_monomer_length"] = None
        self.data["info"]["name"] = self.data["info"]["name"] + "_cryst"
        self.data["symmetries"] = {"all_paths": None,
                                   "unique_paths": None}
        self.paths["export_folder"] = os.path.join(self.paths["export_folder"], "crystal")

        self.monomers = None
        self.ligands = None
        return self

    def set_params(self,
                   min_monomer_length:int,
                   oligomer_levels:int|list[int],
                   ) -> Self:
        """
        Set parameters for crystal processing.
        :param min_monomer_length: Minimum (inclusive) length of chain to be considered a monomer.
        :param oligomer_levels: Oligomerisation level/s to consider.
        :return: Self.
        """
        self.data["crystal"]["min_monomer_length"] = min_monomer_length
        self.data["crystal"]["oligomer_levels"] = oligomer_levels
        return self

    def process(self, force=False) -> Self:
        """
        Processes the crystal through the main pipeline. Requires set_params() to be run beforehand.
        :return: Self.
        """
        log(1, "Processing crystal ({})".format(self.data["info"]["name"]))
        print(self.data)
        print(self.get_full_id())
        self.pass_down()
        self.export()
        self._identyfy_main_elements()
        self.export()
        self._cast_main_elements()
        self.export()
        self._regenerate_crystal()
        self.export()
        self._calculate_oigomerisation_paths()
        self.export()
        self._find_oligomers()
        self.export()
        return self

    def plot(self, paths=False, show=True):
        """
        Plots the Centres of Mass of the monomers and their symmetry mates in a crystal using Matplotlib in
        fractional space.
        """
        fig, ax = fig3D(self, preset="crystal-frac")
        ax.set_title('Crystal {}'.format(self.data["info"]["name"]))

        for n, monomer in enumerate(self.monomers):
            #print(monomer.data["symmetry"]["CoM-frac"])
            col = pymol_colours[n%len(pymol_colours)]
            ax.scatter(*monomer.data["symmetry"]["CoM-frac"], color=col)
            ax.text(*monomer.data["symmetry"]["CoM-frac"], monomer.id, c=col)
            for sym_mon in monomer.sym_elements:
                if sym_mon.data["symmetry"]["operation_n"] == 1:
                    continue

                #print("  -", sym_mon.data["info"]["name"])
                for position in sym_mon.data["symmetry"]["positions"]:
                    displaced_monomer = generate_displaced_copy(
                        monomer.copy(),
                        distance=position,
                        key=sym_mon.data["crystal"]["group_key"],
                        op_n=sym_mon.data["symmetry"]["operation_n"])
                    com = find_com(displaced_monomer)
                    #print("    >", com, position)
                    #cord = [c + p for c,p in zip(cord, position)]
                    ax.scatter(*com, facecolors='none', edgecolors=col)
                    ax.text(*com, sym_mon.data["symmetry"]["operation_n"], c=col)
        if show:
            fig.show()
            input("Press Enter to continue...")
        return fig, ax

    def _identyfy_main_elements(self) -> list[list]:
        """
        Separates chains in model into monomers and ligands, according to the min_monomer_length parameter.
        :return: List of monomers, list of ligands.
        """
        if self.data["crystal"]["min_monomer_length"] is None:
            log("error", "Crystal: missing param: min_monomer_length", raise_exception=True)
        log(2, "Identifying elements ({})".format(self.data["info"]["name"]))
        self.monomers = []
        self.ligands = []
        print(self.get_full_id(), self.data["info"]["name"])
        for chain in self.get_chains():
            c_len = len(chain)
            log(3, chain, c_len, chain.get_full_id(), chain.data["info"]["name"])
            if c_len <= self.data["crystal"]["min_monomer_length"]:
                self.ligands.append(chain)
            else:
                self.monomers.append(chain)
        return [self.monomers, self.ligands]


    def _cast_main_elements(self) -> list[list]:
        """
        Casts monomers and ligands (bi.Chain objects) to their respective classes.
        :return: List of monomers, list of ligands.
        """
        log(1, "Casting main elements ({})".format(self.data["info"]["name"]))
        from .elements import Ligand, Monomer
        for n, mon in enumerate(self.monomers):
            print(mon.data["info"]["name"])
            m = Monomer.cast(mon)
            print(m.data["info"]["name"])
            self.monomers[n] = m
        for n, lig in enumerate(self.ligands):
            l = Ligand.cast(lig)
            self.ligands[n] = l
        log(2, "Monomers: {}".format(self.monomers))
        log(2, "Ligands: {}".format([l.data["info"]["name"] for l in self.ligands]))

        return [self.monomers, self.ligands]


    def _regenerate_crystal(self) -> Self:
        """
        Regenerates the crystal from the given monomers and ligands, and optionally calculates contacts between
        monomers.
        :return:
        """
        log(1, "Regenerating crystal ({})".format(self.data["info"]["name"]))
        from ..visualisation.pymol import PymolScript
        script = PymolScript(name="symmetry_crystal_{}".format(self.data["info"]["name"]),
                             folder=os.path.join(self.paths["export_folder"], "pymol"))
        script.load(self.paths["original"], "original")
        key = self.data["crystal"]["group_key"]
        operations = dictio_space_groups[key]["symops"]
        params = self.data["params"]
        log(2, "Operations:")
        [log(3,o, ">", operations[o]) for o in operations]
        sym_monomers = [] # Fractional
        sym_ligands = [] # Fractional
        log(2, "Monomers ({})".format(len(self.monomers)))
        for monomer in self.monomers:
            log("debug", "Monomer: {}".format(monomer.data["info"]["name"]))
            sym_monomers.extend(monomer.generate_symmetries(self, contacts=True))
        [script.load_entity(entity_to_orth(m.copy(), params)) for m in sym_monomers]

        log(2, "Ligands ({})".format(len(self.ligands)))
        for ligand in self.ligands:
            sym_ligands.extend(ligand.generate_symmetries(self, contacts=False))
        [script.load_entity(entity_to_orth(l.copy(), params)) for l in sym_ligands]

        script.write_script()
        return self


    def _calculate_oigomerisation_paths(self, show=False) -> Self:
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
            self.data["symmetries"]["all_paths"].update(m.data["contacts"]["paths"])
        self.export()
        if show:
            fig.show()
            input("Press Enter to continue...")


        return self


    def _find_oligomers(self, plot=False) -> Self:

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

                    if com2 in point_list:
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


    @staticmethod
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
        from .oligomer import Oligomer
        for path in self.data["symmetries"]["relevant_paths"]:
            print(path)





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

    def add(self, path, contact, reverse=False):
        self.length += 1
        self.data["path"].append({
            "key": path,
            "reverse":reverse
        })
        self.contacts.append(contact)

    def __repr__(self):
        return "<Path in {}: {}, level: {}, closed: {}>".format(
            self.data["name"],
            "-".join([k["key"] for k in self.data["path"]]),
            self.data["o_level"],
            self.data["complete"])


