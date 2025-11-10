
import Bio.PDB as bp
import sys, json
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


    def _calculate_oigomerisation_paths(self) -> Self:
        log(1, "Calculating oligomerisation paths ({})".format(self.data["info"]["name"]))

        path_list = {}
        fig, ax = self.plot(show=False)

        for n, monomer in enumerate(self.monomers):
            log(2, "Monomer: {}".format(monomer.data["info"]["name"]))
            for c in monomer.data["contacts"]["all"]:
                log(3, c["name"])
                com1 = c["monomer1"]["CoM-frac"]
                com2 = c["monomer2"]["CoM-frac"]
                if com2 is None:
                    if c["monomer2"]["operation"] != None:
                        target_mon = [m for m in self.monomers if m.id == c["monomer2"]["id"]][0]
                        print_all_coords(target_mon)
                        print(c["positions"])
                        for position in c["positions"]:
                            displaced_monomer = generate_displaced_copy(
                                target_mon.copy(),
                                distance=position,
                                key=self.data["crystal"]["group_key"],
                                op_n=c["monomer2"]["operation"])
                            com2 = find_com(displaced_monomer)
                            log(4, com1, com2)
                            ax.add_artist(Arrow3D(*zip(com1, com2), color=pymol_colours[n]))
                    else:
                        log("warning", "Monomer2 hs no operation")
                        com2 = find_com([m for m in self.monomers if m.id == c["monomer2"]["id"]][0])
                        log(4, com1, com2)
                        ax.add_artist(Arrow3D(*zip(com1, com2),color=pymol_colours[n]))

        fig.show()
        input("Press Enter to continue...")


        return self


    def _find_oligomers(self) -> Self:

        for m in self.monomers:
            log("debug", m.data["info"]["name"], m.data["symmetry"]["positions"])
            print(m.sym_elements)
        for m in self.ligands:
            log("debug", m.data["info"]["name"], m.data["symmetry"]["positions"])
            print(m.sym_elements)



        for n in self.data["crystal"]["oligomer_levels"]:
            log("debug" "Oligomer level: {}".format(n))

        return self



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
            "a-a":[],
        }


    def add(self, c:dict) -> dict:
        """
        Adds a new atom-atom contact to the list atom-atom contacts. It must be adict with at least the following
        keywords: "is_contact" and "name", additionally a value for "below_threshold" might be required depending on
        the contact method.
        :param c: Dict with the atom-atom contact information.
        """
        self.data["a-a"].append(c)
        return c

    def check_min_contacts(self) -> bool:
        """
        Checks whether this element-element contact meets the required minimum contacts below a certain threshold.
        Requires the min_contacts kwarg on initialisation.
        :return: True if the requirement is met, False otherwise.
        """
        self.data["is_contact"] = False
        pos_dict = {}
        for a in self.data["a-a"]:
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
            self.data["is_contact"], len(self.data["a-a"]),
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













