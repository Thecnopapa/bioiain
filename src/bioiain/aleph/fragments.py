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

    @classmethod
    def from_file(cls, *args, **kwargs):
        self = super().from_file(*args, **kwargs)
        # if not self.has_flag("fragmented", True):
        #     self.fragment(in_place=True)
        self.recover_cvmatrix()
        return self

    def recover_cvmatrix(self):
        pass

    def _fragment_with_aleph(self, force=False, export=False, **kwargs):
        log(2, "Fragmenting structure with ALEPH...")
        if self.has_flag("fragmented", True) and not force:
            log(2, "Fragments already generated!")
            return self._fragments

        if (self.path() is not None) and self.has_flag("fragmented") and not force:
            if os.path.exists(self.path()):
                log(2, "Recovering previously fragmented file...")
                return self.from_file(self.path(), export_folder=self.paths["export_folder"])

        from .core.ALEPH import annotate_pdb_model_with_aleph
        target_path = self.export(minimal=True, target_folder=os.path.join(TEMP_FOLDER, "trash"))
        if target_path is None:
            raise ALEPHError("Minimal cif could not be generated")
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

        log(1, "Fragmented:", target_path)



        fragments = []
        fragmented_ids = graph.vs
        atoms = []
        for n, fraglist in enumerate(fragmented_ids):
            single_atom_ratio = [0, 0]
            cvs = fraglist["reslist"]
            reslist = self.residues()
            target_res= []
            for fres in cvs:
                try:
                    _, model, chain, full_id, resname = fres
                except:
                    #log("warning", fres)
                    _, model, chain, full_id = fres
                    resname = None
                for res in reslist.copy():
                    if res.complex == chain and res.resnum == full_id[1]:
                        if full_id[2].strip() != "":
                            continue # Not implemented
                            print(full_id, res.ca.ins_code)
                            if full_id[2].strip() != res.ca.ins_code:
                                continue
                        if res.resname != resname:
                            if resname is not None:
                                print(res.resseq, res.resnum, res.complex, chain, full_id, res.entity)
                                raise Exception(f"res.resname({res.resname}) != resname({resname})\nfres:{fres}\nres: {res}")
                            else:
                                #log("warning", f"res.resname({res.resname}) != resname({resname})\nfres:{fres}\nres: {res}")
                                pass
                        if len(res.atoms) == 1:
                            single_atom_ratio[0] += 1
                        single_atom_ratio[1] += 1

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
                if single_atom_ratio[0] / single_atom_ratio[1] > 0.5:
                    log("Warning",
                        f"Detected fragment with too many missing side chains ({single_atom_ratio[0] / single_atom_ratio[1]:3.1f}%)")
                    fragment.set_flag("missing_side_chains", True)
                    self.set_flag("missing_side_chains", True)
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
                    log(2, f"fragment {n}: {fragment}", end="\r")

                    fragments.append(fragment)
                if len(fragments) > 0:
                    self._fragments = fragments

        if self._fragments is None or force:
            log(2, "Recalculating fragments...")
            self._fragment_with_aleph()


        return self._fragments



    def _map_cvectors(self, with_ligands=True, vc_mode=None):
        log(1, "Generating CVMatrix for:", self.name(), f"({vc_mode})")
        from . import CVMatrix
        matrix = CVMatrix(self.cvectors(vc_mode=vc_mode), vc_mode=vc_mode, entity=self)
        try:
            matrix.calculate_neighbours()
        except NoNeighboursFound as e:
            log("warning", e)
            self.set_flag("cvmatrix_error", True)
            self._cvmatrix = None
            return None
        matrix.save_fig(attribute="d", save_folder=os.path.join(self.folder(), "cvmaps"))
        if with_ligands:
            matrix.map_ligands(self)
            matrix.save_fig(attribute="dlig", save_folder=os.path.join(self.folder(), "cvmaps"))

        self._cvmatrix = matrix
        self._matrix_vc_mode = vc_mode
        return self._cvmatrix

    def cvmatrix(self, vc_mode=None, **kwargs):
        if self._cvmatrix is None or vc_mode != getattr(self, "_matrix_vc_mode", None):
            self._map_cvectors(vc_mode=vc_mode, **kwargs)
        else:
            log(1, f"Using previously saved CVMatrix ({vc_mode})...")
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

    def show_cvectors(self, execute=True, script=None, vc_mode=None):

        script = self.show_fragments(execute=False, script=script)
        m = self.cvmatrix(vc_mode=vc_mode)
        if m is None:
            log("error", f"No CVMatrix for {self.name()}")
            return None
        for cv in self.cvectors(vc_mode=vc_mode):
            if cv.is_gap:
                script.line("gap_cvectors", coord1=cv.start.coord, coord2=cv.end.coord)
            else:
                script.line("cvectors", coord1=cv.start.coord, coord2=cv.end.coord)

            if m.vc_mode is not None:
                script.line("ca_vc", coord1=cv.res2.ca.coord, coord2=cv.vc.coord)

            script.line("ca_cv", coord1=cv.res2.ca.coord, coord2=cv.start.coord)

            #if cv.fragment is None:
            #    continue
            if cv.closest_opn == 1 or cv.closest_opn is None:
                script.line("neighbours", coord1=cv.vc.coord, coord2=cv.closest.vc.coord)
            else:
                symstart = cv.closest.start.copy().to_frac(self.params()).symop(self.symops(cv.closest_opn), self.params(), position=cv.closest_pos).to_orth(self.params())
                symend = cv.closest.end.copy().to_frac(self.params()).symop(self.symops(cv.closest_opn), self.params(), position=cv.closest_pos).to_orth(self.params())
                script.line("neighbours", coord1=cv.vc.coord, coord2=symstart)
                script.line("cvectors_sym", coord1=symstart, coord2=symend)

        script.color("neighbours", "orange")
        script.color("cvectors", "cyan")
        script.color("gap_cvectors", "purple")
        script.color("cvectors_sym", "green")
        script.color("ca_vc", "yellow")
        script.color("ca_cv", "red")
        print(self.export())
        script.load(self.path(), self.name())
        script.color(self.name(), "white")
        script.hide(self.name())
        script.show(self.name(), "ribbon")
        script.write_script()
        if execute:
            script.execute()
        return script









