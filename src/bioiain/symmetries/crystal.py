
import os, time


from typing_extensions import Self
from copy import deepcopy


from .operations import *
from ..utilities.logging import log
from ..utilities.maths import *
from ..biopython import Model
from .elements import Monomer, Ligand

from ..visualisation import fig3D, pymol_colours, Arrow3D

from .parsing import MissingCrystalError, SuspiciousCrystalError
class Crystal(Model):
    def _init(self, *args, **kwargs) -> Self:
        """
        Initialize the crystal. Executed on casting.
        :param args:
        :param kwargs:
        :return: Self.
        """
        self.force = False
        super()._init(*args, **kwargs)

        if "crystal" not in self.data:
            self.data["crystal"] = {}
        self.data["crystal"]["oligomer_levels"] = None
        self.data["crystal"]["min_monomer_length"] = None
        self.data["info"]["name"] = self.data["info"]["name"] + "_cryst"
        self.data["symmetries"] = {"all_paths": None,
                                   "unique_paths": None}
        self.paths["export_folder"] = os.path.join(self.paths["export_folder"], "crystal")

        self.data["monomers"] = None
        self.data["ligands"] = None
        return self

    def set_crystal_params(self,
                   min_monomer_length:int,
                   min_contacts:int=10,
                   contact_threshold:float|int=15,
                   ) -> Self:
        """
        Set parameters for crystal processing.
        :param min_monomer_length: Minimum (inclusive) length of chain to be considered a monomer.
        :return: Self.
        """
        self.data["crystal"]["min_monomer_length"] = min_monomer_length
        self.data["crystal"]["min_contacts"] = min_contacts
        self.data["crystal"]["contact_threshold"] = contact_threshold
        return self


    def process(self, force=False) -> Self:
        """
        Processes the crystal through the main pipeline. Requires set_params() to be run beforehand.
        :return: Self.
        """
        self.force = force
        log(1, "Processing crystal ({}), FORCE:{}".format(self.data["info"]["name"], self.force))
        try:
            if self.force:
                self.pass_down()
                self.export()
                self._identyfy_main_elements()
                self.export()
                self._regenerate_crystal()
                self.export()

            else:
                self._recover()
                print(self.data["crystal"])
            self.export()
            return self
        except MissingCrystalError:
            return None


    def plot(self, paths=False, show=True):
        """
        Plots the Centres of Mass of the monomers and their symmetry mates in a crystal using Matplotlib in
        fractional space.
        """
        fig, ax = fig3D(self, preset="crystal-frac")
        ax.set_title('Crystal {}'.format(self.data["info"]["name"]))

        for n, monomer in enumerate(self.data["monomers"]):
            monomer = self._restore_monomer(monomer)
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

    def _identyfy_main_elements(self) -> tuple[list, list]:
        """
        Separates chains in model into monomers and ligands, according to the min_monomer_length parameter.
        :return: List of monomers, list of ligands.
        """
        if self.data["crystal"]["min_monomer_length"] is None:
            log("error", "Crystal: missing param: min_monomer_length", raise_exception=True)
        log(2, "Identifying elements ({})".format(self.data["info"]["name"]))

        monomers = []
        ligands = []
        print(self.get_full_id(), self.data["info"]["name"])
        for chain in self.get_chains():
            c_len = len(chain)
            log(3, chain, c_len, chain.get_full_id(), chain.data["info"]["name"])
            if c_len <= self.data["crystal"]["min_monomer_length"]:
                ligands.append(chain)
            else:
                monomers.append(chain)
        self.data["monomers"], self.data["ligands"] = self._cast_main_elements(monomers, ligands)

        return self.data["monomers"], self.data["ligands"]


    def _cast_main_elements(self, monomers, ligands) -> tuple[list, list]:
        """
        Casts monomers and ligands (bi.Chain objects) to their respective classes.
        :return: List of monomers, list of ligands.
        """
        log(1, "Casting main elements ({})".format(self.data["info"]["name"]))
        from .elements import Ligand, Monomer
        mon_ids = []
        lig_ids = []
        self.paths["monomer_folder"] = None
        self.paths["ligand_folder"] = None
        for n, mon in enumerate(monomers):
            m = Monomer.cast(mon)
            m.export()
            mon_ids.append(m.name())
            if n == 0:
                self.paths["monomer_folder"] = m.paths["export_folder"]

        for n, lig in enumerate(ligands):
            l = Ligand.cast(lig)
            l.export()
            lig_ids.append(l.name())
            if n == 0:
                self.paths["ligand_folder"] = l.paths["export_folder"]
        log(2, f"Monomers: {mon_ids}")
        log(2, f"Ligands: {lig_ids}")

        return mon_ids, lig_ids

    def _restore_monomer(self, name):
        return Monomer.recover(name, data_path=os.path.join(self.paths["monomer_folder"], name),
                                  load_structure=True)
    def _restore_ligand(self, name):
        return Ligand.recover(name, data_path=os.path.join(self.paths["ligand_folder"], name),
                                  load_structure=True)

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

        try:
            key = self.data["crystal"]["group_key"]
            operations = dictio_space_groups[key]["symops"]
            params = self.data["params"]
            log(2, "Operations:")
            [log(3, o, ">", operations[o]) for o in operations]
        except KeyError as e:
            raise MissingCrystalError(self)

        sym_monomers = [] # Fractional
        sym_ligands = [] # Fractional
        log(2, "Monomers ({})".format(len(self.data["monomers"])))
        monomers = [self._restore_monomer(m) for m in self.data["monomers"]]
        ligands = [self._restore_ligand(l) for l in self.data["ligands"]]
        for monomer in monomers:
            log("debug", "Monomer: {}".format(monomer.data["info"]["name"]))
            sym_monomers.extend(monomer.generate_symmetries(self, monomers, ligands,
                                                            threshold=self.data["crystal"]["contact_threshold"],
                                                            min_contacts=self.data["crystal"]["min_contacts"],
                                                            contacts=True))
        [script.load_entity(entity_to_orth(m.copy(), params)) for m in sym_monomers]

        log(2, "Ligands ({})".format(len(self.data["ligands"])))
        for ligand in ligands:
            sym_ligands.extend(ligand.generate_symmetries(self, monomers, ligands,
                                                          threshold=self.data["crystal"]["contact_threshold"],
                                                          min_contacts=self.data["crystal"]["min_contacts"],
                                                          contacts=False))
        [script.load_entity(entity_to_orth(l.copy(), params)) for l in sym_ligands]

        script.write_script()
        return self

    def get_oligomers(self, oligomer_levels:int|list[int]):
        self.data["crystal"]["oligomer_levels"] = oligomer_levels
        self._build_oligomers()
        self.export()
        return self



