from .space_groups import dictio_space_groups
import Bio.PDB as bp
import sys
from typing_extensions import Self


from .operations import *
from ..utilities import find_com, print_all_coords
from ..utilities.logging import log
from ..utilities import *
from ..biopython import Chain, Model
from ..visualisation.pymol import quick_display, PymolScript




class CrystalElement(Chain):

    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)
        self.data["symmetry"] = {
            "is_symmetry": False,
            "positions": [],
        }
        self.sym_elements = []



    def generate_symmetries(self, operations, params, key):
        print(self.data["info"]["name"])

        frac_element = entity_to_frac(self.copy(), params)
        frac_element.data = frac_element.data.copy()
        frac_element_com = find_com(frac_element.get_atoms())

        #TODO: Something is not popying right
        for n, operation in operations.items():
            log("debug", "Operation: {}".format(n))
            displaced_element = generate_displaced_copy(frac_element.copy(), key=key, op_n=n).copy()
            displaced_element.data = displaced_element.data.copy()
            displaced_element.data["symmetry"] = {
                "operation_n": n,
                "operation": operation,
                "is_symmetry": True,
                "positions": [],
            }
            displaced_element = displaced_element.copy()
            print(self.data["info"]["name"])
            print(self, displaced_element)
            displaced_element.data["info"]["name"] = "{}_op{}".format(self.data["info"]["name"], n)
            print(self.data["info"]["name"])

            print(displaced_element, [round(c) for c in find_com(frac_element.get_atoms())])
            print(self, [round(c) for c in find_com(self.get_atoms())])

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
                    displaced_element.data["symmetry"]["positions"].append("".join([str(p) for p in position]))

                    if any([p >= 10 for p in position]):
                        print(atom.get_full_id())
                        print("Position:", position)
                        print("Original:", atom.coord)
                        print("Deltas:", deltaX, deltaY, deltaZ)
                        print("New coord:", new_coord)
                        log("error", "Orthogonal position in fractional operation", raise_exception=True)

                    atom.coord = new_coord
                    atom.position = position



            displaced_element.data["symmetry"]["positions"] = list(set(displaced_element.data["symmetry"]["positions"]))
            self.data["symmetry"]["positions"].extend(displaced_element.data["symmetry"]["positions"])
            print(displaced_element.data["symmetry"]["positions"])
            self.sym_elements.append(displaced_element)
            displaced_element.export_data("./exports", displaced_element.data["info"]["name"])
        print(self.data["symmetry"]["positions"])
        print(set(self.data["symmetry"]["positions"]))
        self.data["symmetry"]["positions"] = list(set(self.data["symmetry"]["positions"]))
        print(self.data["symmetry"]["positions"])
        print(self.sym_elements)
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
        self.data["min_monomer_length"] = None
        self.data["oligomer_levels"] = None

        self.monomers = None
        self.ligands = None
        return self

    def set_params(self, data:dict, min_monomer_length, oligomer_levels:int|list[int]) -> Self:
        """
        Set parameters for crystal processing.
        :param data: Data dictionary from parent structure.
        :param min_monomer_length: Minimum (inclusive) length of chain to be considered a monomer.
        :param oligomer_levels: Oligomerisation level/s to consider.
        :return: Self.
        """
        self.data["params"] = data["params"]
        self.data["crystal"] = data["crystal"]
        self.data["info"] = data["info"].copy()
        self.data["min_monomer_length"] = min_monomer_length
        self.data["oligomer_levels"] = oligomer_levels
        return self

    def process(self) -> Self:
        """
        Processes the crystal through the main pipeline. Requires set_params to be run beforehand.
        :return: Self.
        """
        self._identyfy_main_elements()
        self._cast_main_elements()
        self._regenerate_crystal()
        return self

    def _identyfy_main_elements(self) -> list[list[CrystalElement]]:
        """
        Separates chains in model into monomers and ligands, according to the min_monomer_length parameter.
        :return: List of monomers, list of ligands.
        """
        if self.data["min_monomer_length"] is None:
            log("error", "Crystal: missing param: min_monomer_length", raise_exception=True)
        self.monomers = []
        self.ligands = []
        for chain in self.get_chains():
            c_len = len(chain)
            print(chain, c_len)
            if c_len <= self.data["min_monomer_length"]:
                self.ligands.append(chain)
            else:
                self.monomers.append(chain)
        return [self.monomers, self.ligands]


    def _cast_main_elements(self) -> list[list[CrystalElement]]:
        """
        Casts monomers and ligands (bi.Chain objects) to their respective classes.
        :return: List of monomers, list of ligands.
        """
        for n, mon in enumerate(self.monomers):
            m = Monomer.cast(mon)
            self.monomers[n] = m
        for n, lig in enumerate(self.ligands):
            l = Ligand.cast(lig)
            self.ligands[n] = l
        log("debug", "Monomers: {}, Ligands: {}".format(self.monomers, self.ligands))
        return [self.monomers, self.ligands]




    def _regenerate_crystal(self):
        script = PymolScript("symmetry_generation")
        key = self.data["crystal"]["group_key"]
        operations = dictio_space_groups[key]["symops"]
        params = self.data["params"]
        print(operations)
        sym_monomers = [] # Fractional
        sym_ligands = [] # Fractional

        for monomer in self.monomers:
            log("debug", "Monomer: {}".format(monomer))
            sym_monomers.extend(monomer.generate_symmetries(operations, params, key))
        print(sym_monomers)
        [script.load_entity(entity_to_orth(m, params)) for m in sym_monomers]

        for ligand in self.ligands:
            log("debug", "Ligand: {}".format(ligand))
            sym_ligands.extend(ligand.generate_symmetries(operations, params, key))
        [script.load_entity(entity_to_orth(l, params)) for l in sym_ligands]


        for m in sym_monomers:
            log("header", m.data["info"]["name"], m.data["symmetry"]["positions"])


        script.write_script("./exports")
        #script.execute()







    def _find_oligomers(self):

        for n in self.data["oligomer_levels"]:
            log("debug" "Oligomer level: {}".format(n))





class Monomer(CrystalElement):
    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)



class Ligand(CrystalElement):
    def init(self, *args, **kwargs):
        super().init(*args, **kwargs)













