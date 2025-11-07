
import Bio.PDB as bp
import sys
from typing_extensions import Self
from copy import deepcopy

from .space_groups import dictio_space_groups
from .operations import *
from ..utilities import find_com, print_all_coords
from ..utilities.logging import log
from ..utilities import *
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
            "positions": [],
            "CoM-frac": None,
            "CoM-orth": None,
        }
        self.sym_elements = []

    def __repr__(self):
        return "<bi.{} id={} op:{}>".format(
            self.__class__.__name__, self.id, self.data["symmetry"]["operation_n"])

    def __str__(self):
        return "<bi.{} id={} op:{}>".format(
            self.__class__.__name__, self.id, self.data["symmetry"]["operation_n"])



    def generate_symmetries(self, operations, params, key):
        log(3, "Generating symmetries ({})".format(self.data["info"]["name"]))
        log(4, self, "CoM:", [round(c) for c in find_com(self.get_atoms())])
        frac_element = entity_to_frac(self, params)
        frac_element_com = find_com(frac_element.get_atoms())
        self.data["symmetry"]["CoM-frac"] = frac_element_com
        self.data["symmetry"]["CoM-orth"] = find_com(self.get_atoms())

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
            }

            displaced_element.data["info"]["name"] = "{}_op{}".format(self.data["info"]["name"], n)



            for atoms in displaced_element.get_atoms():
                # print(atom, atom.get_full_id(), atom.is_disordered()>0)

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



            displaced_element.data["symmetry"]["positions"] = list(set(
                displaced_element.data["symmetry"]["positions"]))
            displaced_element.data["symmetry"]["CoM-frac"] = find_com(displaced_element.get_atoms())


            self.data["symmetry"]["positions"].extend(displaced_element.data["symmetry"]["positions"])

            self.sym_elements.append(displaced_element)
            log(5, "Element: {}".format(displaced_element.data["info"]["name"]))
            log(5, "Positions:", displaced_element.data["symmetry"]["positions"])

        self.data["symmetry"]["positions"] = list(set(self.data["symmetry"]["positions"]))
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
            sym_monomers.extend(monomer.generate_symmetries(operations, params, key))
        [script.load_entity(entity_to_orth(m.copy(), params)) for m in sym_monomers]

        log(2, "Ligands ({})".format(len(self.ligands)))
        for ligand in self.ligands:
            sym_ligands.extend(ligand.generate_symmetries(operations, params, key))
        [script.load_entity(entity_to_orth(l.copy(), params)) for l in sym_ligands]




        script.write_script()
        #script.execute()


    def _calculate_oigomerisation_paths(self):
        log(1, "Calculating oligomerisation paths ({})".format(self.data["info"]["name"]))

        path_list = {}

        for n, monomer in enumerate(self.monomers):
            log(2, "Monomer: {}".format(monomer.data["info"]["name"]))
            monomer_list = path_list[monomer.id] = {}
            for sym_mon in monomer.sym_elements:
                log(3, "Symmetry: {}".format(sym_mon.data["info"]["name"]))
                for position in sym_mon.data["symmetry"]["positions"]:
                    displaced_monomer = generate_displaced_copy(
                        entity_to_frac(monomer.copy(), monomer.data["params"]),
                        distance=position,
                        key=sym_mon.data["crystal"]["group_key"],
                        op_n=sym_mon.data["symmetry"]["operation_n"])
                    displaced_monomer.data["symmetry"]["position"] = position
                    displaced_monomer.data["info"]["name"] = sym_mon.data["info"]["name"]+"_({})".format(position)
                    log(4, "Position: {}".format(displaced_monomer.data["info"]["name"]))
                    monomer_list[len(monomer_list)] = MonomerContact(monomer, displaced_monomer, mode="min-contacts", threshold=6, min_contacts=3)



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
    def __init__(self, monomer1, monomer2, mode="min-contacts", threshold=None, **kwargs):

        self.data = {
            "is_contact": None,
            "mode": mode,
            "threshold": threshold,
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
            }
        }

        self._calculate_contact(monomer1, monomer2)


    def _calculate_contact(self, monomer1, monomer2):
        log(5, "Calculating contact for: {} - {}".format(monomer1.data["info"]["name"], monomer2.data["info"]["name"]))
        #print(self.data)
        if self.data["mode"] == "min-contacts":
            assert self.data["threshold"] is not None
            assert self.data["kwargs"].get("min_contacts") is not None










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













