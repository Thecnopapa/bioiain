import os, json
import numpy as np

from ..base import BIStructure, BIChain, BIResidue
from ..utilities import *
from ..utilities.exceptions import *


class Fragment(BIChain):
    child_class = BIResidue
    extension = "fragment"

    def __init__(self, *args, fragment_id=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.paths["sub_folder"] = "fragmented/fragments"
        self.data["info"]["fragment_id"] = fragment_id

    def id(self):
        return (self.data["info"]["fragment_id"], self.data["info"]["chain_id"])


    @classmethod
    def from_atoms(cls, atoms, code=None, chain_id=None, fragment_id=None, **kwargs):
        self = super().from_atoms(atoms, code, chain_id, **kwargs)
        self.paths["sub_folder"] = "fragmented/fragments"
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



    def _fragment_with_aleph(self, force=False, export=False, **kwargs):
        log(2, "Fragmenting structure with ALEPH...")
        if (self._fragments is not None) and not force:
            log(2, "Fragments already generated!")
            return self._fragments

        if (self.path() is not None) and self.has_flag("fragmented") and not force:
            if os.path.exists(self.path()):
                log(2, "Recovering previously fragmented file...")
                return self.from_file(self.path(), export_folder=self.paths["export_folder"])

        from .core.ALEPH import annotate_pdb_model_with_aleph
        target_path = self.export(minimal=True, target_folder="/tmp/bioiain/trash")
        self.data["fragments"]["weight"] = "distance_avg"
        self.data["fragments"]["threshold_ah"] = 0.50
        self.data["fragments"]["threshold_bs"] = 0.30
        self.data["fragments"]["peptide_length"] = 3

        try:
            print("### ALEPH start ###")

            graph, _, _, _, _ = annotate_pdb_model_with_aleph(
                target_path,
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
                    if res.resname == resname and res.complex == chain and res.resnum == full_id[1]:
                        #print(res, fres)
                        target_res.append(res)
                        reslist.remove(res)
                        continue
                del res
            del fres
            
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
                log(2, f"fragment {n}: {fragment}")
                fragments.append(fragment)
            else:
                log("warning", f"fragment {n}: no residues")
        #self._atoms = atoms
        self._fragments = fragments

        for a in self.all_atoms():
            if a.get_misc("fragment", None) is None:
                a.set_misc("fragment", None)
                pass

        self.data["fragments"]["n_fragments"] = len(self._fragments)
        self.set_flag("fragmented", True)
        if export:
            self.export()
        #print(self.all_atoms())

        return self._fragments

    def fragments(self, force=False):
        log(1, "Fragmenting structure...")
        if self._fragments is None and not force:
            if self.has_flag("fragmented", True):
                log(2, f"Recovering fragments from cif...")
                atoms_by_fragment = {}
                for a in self.all_atoms():
                    f = a.get_misc("fragment", None)
                    if f is None or f == ".":
                        continue
                    if int(f) not in atoms_by_fragment:
                        atoms_by_fragment[int(f)] = []
                    atoms_by_fragment[int(f)].append(a)
                fragments = []
                for n, fatoms in atoms_by_fragment.items():
                    chain = list(set([r.chain for r in fatoms]))[0]
                    fragment = Fragment.from_atoms(fatoms, code=self.code(), chain_id=chain, fragment_id=n, parent=self)
                    log(2, f"fragment {n}: {fragment}")

                    fragments.append(fragment)
                if len(fragments) > 0:
                    self._fragments = fragments

        if self._fragments is None or force:
            log(2, "Recalculating fragments...")
            self._fragment_with_aleph()


        return self._fragments



    def _map_cvectors(self):
        from . import CVMatrix
        matrix = CVMatrix(self.cvectors())
        matrix.calculate_neighbours()
        self._cvmatrix = matrix
        return self._cvmatrix

    def cvmatrix(self):
        if self._cvmatrix is None:
            print("MAPPING cvmatrix")
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
            script.line("cvectors", coord1=cv.start.coord, coord2=cv.end.coord)
            if cv.closest_opn == 1 or cv.closest_opn is None:
                script.line("neighbours", coord1=cv.start.coord, coord2=cv.closest.start.coord)
            else:
                script.line("neighbours", coord1=cv.start.coord, coord2=cv.closest.start.at(self.symops(cv.closest_opn), self.params(), self.com())[0])
        script.write_script()
        if execute:
            script.execute()
        return script









