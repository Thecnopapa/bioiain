
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
from ..visualisation import fig3D, mpl_colours, pymol_colours
from ..visualisation.pymol import quick_display, PymolScript




class CrystalElement(Chain):

    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)
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



    def generate_symmetries(self, crystal, operations, params, key, contacts=True, threshold=10, min_contacts=1):
        log(3, "Generating symmetries ({})".format(self.data["info"]["name"]))
        log(4, self, "CoM:", [round(c) for c in find_com(self.get_atoms())])
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





class Crystal(Model):
    def init(self, *args, **kwargs) -> Self:
        """
        Initialize the crystal. Executed on casting.
        :param args:
        :param kwargs:
        :return: Self.
        """
        super().init(*args, **kwargs)
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
        :param data: Data dictionary from parent structure.
        :param min_monomer_length: Minimum (inclusive) length of chain to be considered a monomer.
        :param oligomer_levels: Oligomerisation level/s to consider.
        :return: Self.
        """
        self.data["crystal"]["min_monomer_length"] = min_monomer_length
        self.data["crystal"]["oligomer_levels"] = oligomer_levels
        self.paths["crystal_folder"] = crystal_folder
        return self

    def process(self) -> Self:
        """
        Processes the crystal through the main pipeline. Requires set_params to be run beforehand.
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

    def plot(self):
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

        fig.show()
        input("Press Enter to continue...")

    def _identyfy_main_elements(self) -> list[list[CrystalElement]]:
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


    def _cast_main_elements(self) -> list[list[CrystalElement]]:
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


    def _regenerate_crystal(self):
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
            sym_monomers.extend(monomer.generate_symmetries(self, operations, params, key, contacts=True))
        [script.load_entity(entity_to_orth(m.copy(), params)) for m in sym_monomers]

        log(2, "Ligands ({})".format(len(self.ligands)))
        for ligand in self.ligands:
            sym_ligands.extend(ligand.generate_symmetries(self, operations, params, key, contacts=False))
        [script.load_entity(entity_to_orth(l.copy(), params)) for l in sym_ligands]




        script.write_script()
        #script.execute()


    def _calculate_oigomerisation_paths(self):
        log(1, "Calculating oligomerisation paths ({})".format(self.data["info"]["name"]))

        path_list = {}

        for n, monomer in enumerate(self.monomers):
            log(2, "Monomer: {}".format(monomer.data["info"]["name"]))
            [log(3, c["name"]) for c in monomer.data["contacts"]["all"]]



    def _find_oligomers(self):

        for m in self.monomers:
            log("debug", m.data["info"]["name"], m.data["symmetry"]["positions"])
            print(m.sym_elements)
        for m in self.ligands:
            log("debug", m.data["info"]["name"], m.data["symmetry"]["positions"])
            print(m.sym_elements)



        for n in self.data["crystal"]["oligomer_levels"]:
            log("debug" "Oligomer level: {}".format(n))



class MonomerContact(object):
    def __init__(self, monomer1, monomer2, mode="min-contacts", threshold=None, min_contacts=None, **kwargs):

        self.data = {
            "name":"Contact: ({}-{}): Unprocessed".format(monomer1.data["info"]["name"], monomer2.data["info"]["name"]),
            "is_contact": None,
            "mode": mode,
            "threshold": threshold,
            "min_contacts": min_contacts,
            "kwargs": kwargs,
            "monomer1": {
                "name": monomer1.data["info"]["name"],
                "id": monomer1.id,
                "operation": monomer1.data["symmetry"]["operation_n"],
                "is_symmetry": monomer1.data["symmetry"]["is_symmetry"],
                "position": monomer1.data["symmetry"].get("position", None),
            },
            "monomer2": {
                "name": monomer2.data["info"]["name"],
                "id": monomer2.id,
                "operation": monomer2.data["symmetry"]["operation_n"],
                "is_symmetry": monomer2.data["symmetry"]["is_symmetry"],
                "position": monomer2.data["symmetry"].get("position", None),
            },
            "a-a":[],
        }

        #self._calculate_contact(monomer1, monomer2)
        #print(json.dumps(self.data, indent=4))

    def add(self, c):
        self.data["a-a"].append(c)

    def check_min_contacts(self):
        count = 0
        for a in self.data["a-a"]:
            if a["below_threshold"]:
                count += 1
                if count >= self.data["min_contacts"]:
                    self.data["is_contact"] = True
                    self.data["name"] = repr(self)
                    return True
        self.data["is_contact"] = False
        self.data["name"] = repr(self)
        return False


    def __repr__(self):
        return "Contact: ({} <-> {}): {}, N:{}, T:{}, min:{}".format(
            self.data["monomer1"]["name"], self.data["monomer2"]["name"],
            self.data["is_contact"], len(self.data["a-a"]), self.data["threshold"], self.data["min_contacts"])

    def _calculate_contact(self, monomer1, monomer2):
        log(5, "Calculating contact for: {} - {}".format(monomer1.data["info"]["name"], monomer2.data["info"]["name"]))


        contacts = []

        if self.data["mode"] == "min-contacts":
            assert self.data["threshold"] is not None
            assert self.data["kwargs"].get("min_contacts") is not None

            print_all_coords(monomer1)
            print_all_coords(monomer2)

            monomer1 = entity_to_orth(monomer1.copy(), monomer1.data["params"])
            monomer2 = entity_to_orth(monomer2.copy(), monomer1.data["params"])


            for a1 in monomer1.get_atoms():
                for a2 in monomer2.get_atoms():
                    d = d2(a1.coord, a2.coord, root=False)
                    print(a1, a2, d, end="\r")
                    if d < self.data["threshold"]**2:

                        contacts.append({
                            "atom1": a1.get_full_id(),
                            "atom2": a2.get_full_id(),
                            "distance": d,
                            "below_threshold": True,
                            "threshold": self.data["threshold"],
                        })

        self.data["contacts"] = contacts










class Monomer(CrystalElement):
    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)



class Ligand(CrystalElement):
    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)













