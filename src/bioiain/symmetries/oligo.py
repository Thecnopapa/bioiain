from .space_groups import dictio_space_groups
import Bio.PDB as bp
import sys

from .operations import *
from ..utilities import find_com, print_all_coords
from ..utilities.logging import log
from ..utilities import *
from ..biopython import Chain, Model
from ..visualisation.pymol import quick_display, PymolScript




class Crystal(Model):
    def init(self, *args, **kwargs):
        self.data["min_monomer_length"] = None
        self.data["oligomer_levels"] = None

        self.monomers = None
        self.ligands = None

    def set_params(self, data:dict, min_monomer_length, oligomer_levels):
        self.data["params"] = data["params"]
        self.data["crystal"] = data["crystal"]
        self.data["min_monomer_length"] = min_monomer_length
        self.data["oligomer_levels"] = oligomer_levels

    def process(self):
        self._identyfy_main_elements()
        self._cast_main_elements()
        self._regenerate_crystal()

    def _identyfy_main_elements(self):
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


    def _cast_main_elements(self):
        for n, mon in enumerate(self.monomers):
            m = Monomer.cast(mon)
            self.monomers[n] = m
        for n, lig in enumerate(self.ligands):
            l = Ligand.cast(lig)
            self.ligands[n] = l
        log("debug", "Monomers: {}, Ligands: {}".format(self.monomers, self.ligands))


    def _find_oligomers(self):

        for n in self.data["oligomer_levels"]:
            log("debug" "Oligomer level: {}".format(n))




    def _regenerate_crystal(self):
        script = PymolScript("symmetry_generation")
        key = self.data["crystal"]["group_key"]
        operations = dictio_space_groups[key]["symops"]
        params = self.data["params"]
        print(operations)
        sym_monomers = []
        sym_ligands = []


        for monomer in self.monomers:
            log("debug", "Monomer: {}".format(monomer))
            frac_monomer = entity_to_frac(monomer.copy(), params)
            frac_monomer_com =  find_com(frac_monomer.get_atoms())
            for n, operation in operations.items():
                log("debug", "Operation: {}".format(n))
                displaced_monomer = generate_displaced_copy(frac_monomer, distance=99.5, key=key, op_n=n)
                print(displaced_monomer, [round(c) for c in find_com(displaced_monomer.get_atoms())])
                print(monomer, [round(c) for c in find_com(monomer.get_atoms())])

                for atoms in displaced_monomer.get_atoms():
                    #print(atom, atom.get_full_id(), atom.is_disordered()>0)

                    if not atoms.is_disordered()>0:
                        atoms = [atoms]

                    for atom in atoms:

                        deltaX = ((atom.coord[0] - frac_monomer_com[0]) % 1) - 0.5
                        deltaY = ((atom.coord[1] - frac_monomer_com[1]) % 1) - 0.5
                        deltaZ = ((atom.coord[2] - frac_monomer_com[2]) % 1) - 0.5

                        new_coordX = frac_monomer_com[0] + deltaX
                        new_coordY = frac_monomer_com[1] + deltaY
                        new_coordZ = frac_monomer_com[2] + deltaZ
                        new_coord = [new_coordX, new_coordY, new_coordZ]

                        position = [(n_coord - d_coord + 99.5) for n_coord, d_coord in zip(new_coord, atom.coord)]
                        for p in position:
                            assert p % 1 == 0
                        position = tuple([int(p) for p in position])

                        if any([p >= 10 for p in position]):
                            print(atom.get_full_id())
                            print("Position:", position)
                            print("Original:", atom.coord)
                            print("Deltas:", deltaX, deltaY, deltaZ)
                            print("New coord:", new_coord)
                            log("error", "Orthogonal position in fractional operation", raise_exception=True)

                        atom.coord = new_coord
                        atom.position = position

                entity_to_orth(displaced_monomer, params)
                sym_monomers.append(displaced_monomer)
                script.load_entity(displaced_monomer)
        script.write_script("./exports")
        script.execute()


















        #quick_display(self.fractional)
        #quick_display(self)




class CrystalElement(Chain):
    pass


class Monomer(CrystalElement):
    pass



class Ligand(CrystalElement):
    pass














