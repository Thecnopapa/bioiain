
import Bio.PDB as bp
import sys, json

from docutils.utils.math.mathml_elements import mover
from typing_extensions import Self
from copy import deepcopy

from .space_groups import dictio_space_groups
from .operations import *
from ..utilities import find_com, print_all_coords
from ..utilities.logging import log
from ..utilities import *
from ..utilities.maths import *
from ..biopython import Chain, Model
from ..visualisation import fig3D, mpl_colours, pymol_colours, Arrow3D
from ..visualisation.pymol import quick_display, PymolScript







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
        self.paths["crystal_folder"] = None
        self.data["info"]["name"] = self.data["info"]["name"] + "_cryst"
        self.data["symmetries"] = {"all_paths": None,
                                   "unique_paths": None}

        self.monomers = None
        self.ligands = None
        return self

    def set_params(self,
                   min_monomer_length:int,
                   oligomer_levels:int|list[int],
                   crystal_folder:str="./crystals"
                   ) -> Self:
        """
        Set parameters for crystal processing.
        :param min_monomer_length: Minimum (inclusive) length of chain to be considered a monomer.
        :param oligomer_levels: Oligomerisation level/s to consider.
        :param crystal_folder: Folder to store crystal exports.
        :return: Self.
        """
        self.data["crystal"]["min_monomer_length"] = min_monomer_length
        self.data["crystal"]["oligomer_levels"] = oligomer_levels
        self.paths["crystal_folder"] = crystal_folder
        return self

    def process(self) -> Self:
        """
        Processes the crystal through the main pipeline. Requires set_params() to be run beforehand.
        :return: Self.
        """
        log(1, "Processing crystal ({})".format(self.data["info"]["name"]))
        print(self.data)
        print(self.get_full_id())
        self.pass_down()
        self._identyfy_main_elements()
        self._cast_main_elements()
        self._regenerate_crystal()
        self._calculate_oigomerisation_paths()
        self._find_oligomers()
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
        script = PymolScript(name="symmetry_crystal_{}".format(self.data["info"]["name"]),
                             folder=self.paths["crystal_folder"])
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


    def _find_oligomers(self) -> Self:

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


        for n, path in enumerate(unique_paths):
            if path["o_level"] == 2:
                continue
            log(3, "Path: {}, level: {}".format(n, path["o_level"]) )
            log(4, "Steps: {}".format(path["steps"]))
            starting_monomer = path["steps"][0][0]
            log(4, "Starting monomer: {}".format(starting_monomer))


            fig, ax = self.plot(show=False)

            point_list = []
            operation_list = []
            current_pos = [0, 0, 0]

            coms = deepcopy(o_coms)
            print("COMS:", coms)
            point_list.append(coms[starting_monomer])

            for s, step in enumerate(path["path"]):
                step_info = all_paths[step["key"]]
                #print("###")
                #print(step_info)
                print("###")
                print(step)
                #print("###")
                op_n = step_info["monomer2"]["operation"]
                key = self.data["crystal"]["group_key"]
                params = self.data["params"]
                pos = step_info["position"]

                reverse = step["reverse"]

                # if reverse:
                #     new_pos = [c-p for p,c in zip(pos, current_pos)]
                # else:
                #     new_pos = [c+p for p,c in zip(pos, current_pos)]

                #print("new_pos:", new_pos)
                print("pos:", pos)


                # print(dict(
                #     op_n=op_n,
                #     key=key,
                #     params=params,
                #     position=position,
                #     reverse=reverse,
                # ))

                id1 = step_info["monomer1"]["id"]
                id2 = step_info["monomer2"]["id"]
                d = [1/p if p != 0 else 0 for p in pos]
                #if reverse:
                #    d = [-p for p in pos]

                if reverse:
                    id1, id2 = id2, id1

                print(id1, "-->", id2)




                # Maybe log all operations and carry them every time, but should be the same
                print("current pos:", current_pos)
                print(coms[id2], "-->", end=" ")
                current_pos = [c + p for c, p in zip(current_pos, d)]
                print(current_pos)
                dcoms = {k: coord_add(v,d, True) for k, v in coms.items()}
                ax.scatter(*dcoms[id2])
                ax.text()
                coms = {k: coord_operation(v, key, op_n, current_pos) for k, v in dcoms.items()}
                moving_com = [c+p for c, p in zip(coms[id2], d)]
                print(moving_com)



                #moving_com = [p%1 for p in coms[id2]]

                point_list.append(moving_com)
                ax.add_artist(Arrow3D(*zip(point_list[-2], point_list[-1]), color="black"))

                plane = {"x": [0,10], "y": [0,10], "z": [0,10]}
                if pos[0] != 0:
                    plane["x"][0] = int(1 / pos[0] * 10)
                if pos[1] != 0:
                    plane["y"][0] = int(1 / pos[1] * 10)
                if pos[2] != 0:
                    plane["z"][0] = int(1 / pos[2] * 10)



                # Check for early circle closures
                # Remove if so

                # Then classify linear/circular based on final point


            fig.show()
            input("Press Enter to continue...")



            print(point_list)
            print(operation_list)





        self.data["symmetries"]["unique_paths"] = unique_paths


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
                                 if (k.startswith(last_id) or k.endswith(last_id)) and k != last_key}
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




class CrystalElement(Chain):

    def _init(self, *args, **kwargs):
        super()._init(*args, **kwargs)
        self.data["symmetry"] = {
            "operation": None,
            "operation_n": None,
            "is_symmetry": False,
            "positions": None,
            "CoM-frac": None,
            "CoM-orth": None,
            "contacts": None,
        }
        self.data["contacts"] = {
            "all": [],
            "paths": None,
        }
        self.sym_elements = []

    def __repr__(self):
        return "<bi.{} id={} op:{}>".format(
            self.__class__.__name__, self.id, self.data["symmetry"]["operation_n"])

    def __str__(self):
        return "<bi.{} id={} op:{}>".format(
            self.__class__.__name__, self.id, self.data["symmetry"]["operation_n"])



    def generate_symmetries(self, crystal:Crystal,contacts:bool=True,
                            threshold:int|float=10, min_contacts:int=1,
                            contact_method:str="min-contacts") -> list[Self]:
        """
        Generates symmetries of crystal elements from a given symmetry element, and (optionally) calculates contacts
        (theoretically between monomers).
        A distance threshold in Angstroms and a minimum number of contacts is required for contact calculations.
        The symmetry Elements generated contain one structure per operation, which is split across different positions.
        The contacts generated are stored at a monomer level.
        :param crystal: Crystal object.
        :param contacts: Whether to calculate contacts between monomers.
        :param threshold: Threshold in Angstroms to consider a CA-CA contact.
        :param min_contacts: Minimum contacts to consider a monomer-element contact.
        :param contact_method: Method to determine contacts between elements.
        :return: List of generated symmetry Elements
        """
        log(3, "Generating symmetries ({})".format(self.data["info"]["name"]))
        log(4, self, "CoM:", [round(c) for c in find_com(self.get_atoms())])
        params = crystal.data["params"]
        key = crystal.data["crystal"]["group_key"]
        operations = dictio_space_groups[key]["symops"]

        frac_element = entity_to_frac(self, params)
        frac_element_com = find_com(frac_element.get_atoms())
        self.data["symmetry"]["CoM-frac"] = frac_element_com
        self.data["symmetry"]["CoM-orth"] = find_com(self.get_atoms())
        self.data["symmetry"]["positions"] = []
        if contacts:
            self.data["contacts"]["all"] = []
            self.data["contacts"]["threshold"] = threshold
            self.data["contacts"]["min_contacts"] = min_contacts

        for n, operation in operations.items():

            log(4, "Operation: {}".format(n))
            displaced_element = generate_displaced_copy(frac_element.copy(), key=key, op_n=n).copy()
            displaced_element.data = deepcopy(displaced_element.data)
            displaced_element.data["symmetry"] = {
                "operation_n": n,
                "operation": operation,
                "is_symmetry": True,
                "positions": [],
                "CoM-frac": None,
                "CoM-orth": None,
                "contacts": {}
            }
            displaced_element.data.pop("contacts")

            displaced_element.data["info"]["name"] = "{}_op{}".format(self.data["info"]["name"], n)

            if contacts:
                contact = {m.id:MonomerContact(m, displaced_element,
                                               contact_method=contact_method,
                                               threshold=threshold,
                                               min_contacts=min_contacts) for m in crystal.monomers if not (
                    n == 1 and displaced_element.id == m.id
                )}



            for atoms in displaced_element.get_atoms():

                if not atoms.is_disordered() > 0:
                    atoms = [atoms]

                for atom in atoms:

                    deltaX = ((atom.coord[0] - frac_element_com[0]) % 1) - 0.5
                    deltaY = ((atom.coord[1] - frac_element_com[1]) % 1) - 0.5
                    deltaZ = ((atom.coord[2] - frac_element_com[2]) % 1) - 0.5

                    new_coordX = frac_element_com[0] + deltaX
                    new_coordY = frac_element_com[1] + deltaY
                    new_coordZ = frac_element_com[2] + deltaZ
                    new_coord = [new_coordX, new_coordY, new_coordZ]

                    position = [(n_coord - d_coord + 99.5) for n_coord, d_coord in zip(new_coord, atom.coord)]
                    for p in position:
                        assert p % 1 == 0
                    position = tuple([int(p) for p in position])
                    displaced_element.data["symmetry"]["positions"].append(position)

                    if any([p >= 10 for p in position]):
                        print(atom.get_full_id())
                        print("Position:", position)
                        print("Original:", atom.coord)
                        print("Deltas:", deltaX, deltaY, deltaZ)
                        print("New coord:", new_coord)
                        log("error", "Orthogonal position in fractional operation", raise_exception=True)

                    atom.coord = new_coord
                    atom.position = position

                    if contacts and atom.id == "CA":
                        for m in crystal.monomers:
                            if n == 1 and displaced_element.id == m.id:
                                continue
                            for a in m.get_atoms():
                                if not a.id == "CA":
                                    continue
                                d = get_fractional_distance(a.coord, atom.coord, self.data["params"])
                                print(atom.parent.id[1], atom.parent.resname, "\t", a.parent.id[1], a.parent.resname, "\t", d, end="\t")

                                if d <= threshold**2:
                                    print("true", end="\r")
                                    contact[m.id].add({
                                        "atom1": a.get_full_id(),
                                        "atom2": atom.get_full_id(),
                                        "distance": d,
                                        "below_threshold": True,
                                        "threshold": threshold,
                                        "position": atom.position,
                                    })
                                else:
                                    print("false", end="\r")

            if contacts:

                for m in crystal.monomers:
                    if n == 1 and displaced_element.id == m.id:
                        continue
                    displaced_element.data["symmetry"]["contacts"][m.id] = contact[m.id].data
                    if contact[m.id].check_min_contacts():
                        m.data["contacts"]["all"].append(contact[m.id].data)





            displaced_element.data["symmetry"]["positions"] = list(set(
                displaced_element.data["symmetry"]["positions"]))
            displaced_element.data["symmetry"]["CoM-frac"] = find_com(displaced_element.get_atoms())


            self.data["symmetry"]["positions"].extend(displaced_element.data["symmetry"]["positions"])

            self.sym_elements.append(displaced_element)
            log(5, "Element: {}".format(displaced_element.data["info"]["name"]))
            log(5, "Positions:", displaced_element.data["symmetry"]["positions"])

        if contacts:
            log(4, "Elements in contact: {}".format(len(self.data["contacts"]["all"])))
            for c in self.data["contacts"]["all"]:
                log(5, c["name"])

        self.data["symmetry"]["positions"] = list(set(self.data["symmetry"]["positions"]))
        self.export(structure=False)
        return self.sym_elements







class MonomerContact(object):
    """
    Class to calculate and store element-element contact information.
    :param monomer1: First monomer (certain values, such as the parameter dict will be obtained from this entity).
    :param monomer2: Second monomer.
    :param kwargs: Other parameters to store in the contact such as the threshold used.
    """
    def __init__(self, monomer1, monomer2, **kwargs):
        self.data = {
            "name":"Contact: ({}-{}): Unprocessed".format(monomer1.data["info"]["name"], monomer2.data["info"]["name"]),
            "is_contact": None,
            "positions": [],
            "position": None,
            "kwargs": kwargs,
            "monomer1": {
                "name": monomer1.data["info"]["name"],
                "id": monomer1.id,
                "operation": monomer1.data["symmetry"]["operation_n"],
                "is_symmetry": monomer1.data["symmetry"]["is_symmetry"],
                "position": monomer1.data["symmetry"].get("position", None),
                "CoM-frac": monomer1.data["symmetry"]["CoM-frac"],
            },
            "monomer2": {
                "name": monomer2.data["info"]["name"],
                "id": monomer2.id,
                "operation": monomer2.data["symmetry"]["operation_n"],
                "is_symmetry": monomer2.data["symmetry"]["is_symmetry"],
                "position": monomer2.data["symmetry"].get("position", None),
                "CoM-frac": monomer2.data["symmetry"]["CoM-frac"],

            },
        }
        self.a_a = []


    def add(self, c:dict) -> dict:
        """
        Adds a new atom-atom contact to the list atom-atom contacts. It must be adict with at least the following
        keywords: "is_contact" and "name", additionally a value for "below_threshold" might be required depending on
        the contact method.
        :param c: Dict with the atom-atom contact information.
        """
        self.a_a.append(c)
        return c

    def check_min_contacts(self) -> bool:
        """
        Checks whether this element-element contact meets the required minimum contacts below a certain threshold.
        Requires the min_contacts kwarg on initialisation.
        :return: True if the requirement is met, False otherwise.
        """
        self.data["is_contact"] = False
        pos_dict = {}
        for a in self.a_a:
            if a["below_threshold"]:
                if str(["position"]) in pos_dict:
                    pos_dict[str(["position"])] += 1
                else:
                    pos_dict[str(["position"])] = 1
                if pos_dict[str(["position"])]  >= self.data["kwargs"].get("min_contacts", 1):
                    self.data["is_contact"] = True
                    self.data["positions"].append(a["position"])
                    self.data["name"] = repr(self)
        if self.data["is_contact"]:
            self.data["positions"] = list(set(self.data["positions"]))
            return True
        else:
            self.data["name"] = repr(self)
            return False


    def __repr__(self):
        return "Contact: ({} <-> {}): {}, N:{}, T:{}, min:{}".format(
            self.data["monomer1"]["name"], self.data["monomer2"]["name"],
            self.data["is_contact"], len(self.a_a),
            self.data["kwargs"].get("threshold", None), self.data["kwargs"].get("min_contacts", None))









class Monomer(CrystalElement):
    def _init(self, *args, **kwargs):
        super()._init(*args, **kwargs)

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)



class Ligand(CrystalElement):
    def _init(self, *args, **kwargs):
        super()._init(*args, **kwargs)

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)













