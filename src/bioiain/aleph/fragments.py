import os, json
import numpy as np

from . import CVMatrix
from ..base import BIStructure, BIChain, BIResidue
from ..utilities import *
from ..utilities.exceptions import *


class Fragment(BIChain):
    child_class = BIResidue
    extension = "fragment"

    def __init__(self, *args, fragment_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.paths["sub_folder"] = "fragmented/fragments"


    @classmethod
    def from_atoms(cls, atoms, code=None, chain_id=None, fragment_id=None, **kwargs):
        self = super().from_atoms(atoms, code, chain_id, **kwargs)
        self.data["info"]["fragment_id"] = f"F{fragment_id}"
        self.set_name(fragment_id, append=True)
        for a in self.atoms():
            a.set_misc("fragment", fragment_id)
        return self






class FragmentedStructure(BIStructure):
    child_class = Fragment
    extension = "fstructure"

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.paths["sub_folder"] = "fragmented"
        self.data["fragments"] = {}
        self._cvmatrix = None
        self._fragments = None



    def fragment_with_aleph(self, force=False, export=False):
        from .core.ALEPH import annotate_pdb_model_with_aleph

        if (self._fragments is not None) and not force:
            return self._fragments

        self.export(minimal=True)
        self.data["fragments"]["weight"] = "distance_avg"
        self.data["fragments"]["threshold_ah"] = 0.50
        self.data["fragments"]["threshold_bs"] = 0.30
        self.data["fragments"]["peptide_length"] = 3

        try:
            print("### ALEPH start ###")

            graph, _, _, _, _ = annotate_pdb_model_with_aleph(
                self.paths["minimal"],
                weight=self.data["fragments"]["weight"],
                strictness_ah=self.data["fragments"]["threshold_ah"],
                strictness_bs=self.data["fragments"]["threshold_bs"],
                peptide_length=self.data["fragments"]["peptide_length"],
                write_pdb=False,
            )
            print("### ALEPH end ###")
        except Exception as e:
            print("### ALEPH failed ###")
            log("warning", e)
            raise ALEPHError(e)

        log(1, "Fragmented:", self.paths["minimal"],)



        fragments = []
        fragmented_ids = graph.vs
        atoms = []
        for n, fraglist in enumerate(fragmented_ids):
            cvs = fraglist["reslist"]
            reslist = self.residues()
            target_res= []
            for fres in cvs:
                try:
                    _, model, chain, full_id, resname = fres
                except:
                    log("warning", fres)
                    _, model, chain, full_id = fres
                    resname = None
                for res in reslist.copy():
                    if res.resname == resname and res.chain == chain and res.resnum == full_id[1]:
                        #print(res, fres)
                        target_res.append(res)
                        reslist.remove(res)
                        continue
                del res
            del fres
            log(2, f"fragment {n}: {len(target_res)}")
            if len(target_res) == 0:
                log("warning", f"fragment {n}: no residues")
                break

            fatoms = []
            for res in target_res:
                fatoms.extend(res.atoms)
                atoms.extend(res.atoms)
            chain = list(set([r.chain for r in fatoms]))
            if len(chain) > 0:
                chain = chain[0]
                fragment = Fragment.from_atoms(fatoms, code=self.code(), chain_id=chain, fragment_id=n, parent=self)
                if export:
                    fragment.export()
                fragments.append(fragment)
            else:
                log("warning", f"fragment {n}: no residues")
        #self._atoms = atoms
        self._fragments = fragments
        #for a in self.all_atoms():
            #if a.get_misc("fragment", None) is None:
                #self.remove_atom(a)
        self.data["fragments"]["n_fragments"] = len(self._fragments)
        if export:
            self.export()
        #print(self.all_atoms())

        return self._fragments

    def fragments(self):
        if self._fragments is None:
            self.fragment_with_aleph()

        return self._fragments



    def _map_cvectors(self):
        matrix = CVMatrix(self.cvectors())
        matrix.calculate_neighbours()
        self._cvmatrix = matrix
        return self._cvmatrix

    def cvmatrix(self):
        if self._cvmatrix is None:
            self._map_cvectors()
        return self._cvmatrix

    def show(self):
        for a in self.all_atoms():
            a.set_bfactor(a.get_misc("fragment"))
        super().show()

    def show_fragments(self, execute=True, script=None):
        if script is None:
            from ..visualisation.pymol import PymolScript
            script = PymolScript(self.name())
        for fragment in self.fragments():
            for a in fragment.all_atoms():
                a.set_bfactor(a.get_misc("fragment"))
            script.load(fragment.export(), "frag_"+fragment.name())
        script.group("frag_", "fragments")
        script.spectrum("(all)")
        script.orient()
        script.write_script()
        if execute:
            script.execute()
        return script

    def show_cvectors(self, execute=True, script=None):
        script = self.show_fragments(execute=False, script=script)
        self.cvmatrix()
        for cv in self.cvectors():
            script.line("cvectors", coord1=cv.start, coord2=cv.end)
            script.line("neighbours", coord1=cv.start, coord2=cv.closest.start)
        script.write_script()
        if execute:
            script.execute()
        return script









