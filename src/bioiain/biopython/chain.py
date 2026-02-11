import os, sys, json
import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .residue import Residue, DResidue
from ..utilities.sequences import d3to1
from ..utilities.logging import log
import numpy as np
import math


class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    child_class = Residue
    disordered_child_class = DResidue

    def _init(self, *args, **kwargs):
        self.data["sequence"] = None
        self._kdtree = None


    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def residues():
        pass

    def atoms(self, ca_only=False, hetatm=False, force=False, group_by_residue=False, disordered=False, **kwargs):
        from .imports import read_mmcif
        from .atom import BIAtom

        _params = {ca_only:ca_only, hetatm:hetatm, group_by_residue:group_by_residue, disordered:disordered}

        if not force:
        #    if hasattr(self, "_atoms_params"):
        #        for k, v in _params.items():
        #            if self._atoms_params.get(k, None) != v:
        #                force= True
        #                break
        #    else:
        #        force = True
            pass

        if not hasattr(self, "_atoms"):
            force = True
        elif self._atoms is None:
            force = True


        if force:
            print("Reading atoms from CIF")
            atoms = read_mmcif(self.paths["self"], subset=["_atom_site"])("_atom_site")
            atoms = [BIAtom(a) for a in atoms]

            self._atoms = atoms
            self._atoms_params = {
                "ca_only": ca_only,
                "hetatm": hetatm,
                "group_by_residue": group_by_residue,
                "disordered": disordered,
            }

        atoms = self._atoms
        if not disordered:
            atoms = self._fix_disordered(atoms)

        if not hetatm:
            atoms = [a for a in atoms if a.type != "HETATM"]
        if ca_only:
            atoms = [a for a in atoms if a.name == "CA"]
        if group_by_residue:
            atoms_by_res = {}
            residues_with_ca = []
            for atom in atoms:
                if atom.resnum in atoms_by_res:
                    atoms_by_res[atom.resnum].append(atom)
                else:
                    atoms_by_res[atom.resnum] = [atom]
                    
                if atom.name == "CA":
                    residues_with_ca.append(atom.resnum)
            for key in list(atoms_by_res.keys()):
                try:
                    residues_with_ca.remove(key)
                except:
                    atoms_by_res.pop(key)

            atoms = atoms_by_res

        return atoms


    @staticmethod
    def _fix_disordered(atoms):
        fixed_atoms = []
        for atom in atoms:
            if atom.disordered:
                for a in fixed_atoms:
                    if a.id == atom.id:
                        a.doppelgangers.append(atom)
                        a.favourite = True
                        atom.favourite = False
                        atom.doppelgangers = None
                        break
                if atom.favourite:
                    fixed_atoms.append(atom)
            else:
                fixed_atoms.append(atom)
        return fixed_atoms



    def get_sequence(self, force=False):

        atoms = self.atoms(ca_only=True)

        if self.data["sequence"] is not None and not force:
            return self.data["sequence"]
        self.data["sequence"] = ""
        for atom in atoms:
            try:
                #print(atom["label_comp_id"])
                r = d3to1[atom.resname]
                #print(r)
                self.data["sequence"] += r
            except KeyError as e:
                log("debug", f"Residue {e} not recognised")
                #input("press enter")
        return self.data["sequence"]


    def compute_sasa(self, **kwargs):
        print("computing SASA...")
        from .SASA import SASA
        sasa = SASA(**kwargs)
        return sasa.compute(self, **kwargs)

    def get_surface_residues(self, threshold = 50, ball_radius=1.40, force=False, reset_other=True):

        if "surface" in self.data and not force:
            if str(threshold) in self.data["surface"]:
                if self.data["surface"][str(threshold)]["reslist"] is not None:
                    if self.data["surface"][str(threshold)].get("ball_radius", None) == ball_radius:
                        print("Using precalculated SASA")
                        return self.data["surface"][str(threshold)]["reslist"]
                    else:
                        force = True

        atoms_by_res = self.atoms(group_by_residue=True)
        try:
            if list(atoms_by_res.values())[0][0].get_misc("SASA") is not None:
                raise KeyError
        except:
            self.compute_sasa(ball_radius=ball_radius, force=force)

        surface_res_ids = []

        for resn, atom_group in atoms_by_res.items():
            res_asa = sum([float(a.get_misc("SASA")) for a in atom_group]) / len(atom_group)
            #print(resn, res_asa)
            if res_asa >= threshold:
                surface_res_ids.append(resn)

        percetage_outer = len(atoms_by_res) / len(surface_res_ids)

        if "surface" not in self.data or reset_other:
            self.data["surface"] = {str(threshold): {}}

        self.data["surface"][str(threshold)] = {
        "reslist": surface_res_ids,
        "ball_radius": ball_radius,
        "percetage_outer": percetage_outer,
        }

        self.export(include_misc=True)
        return self.data["surface"][str(threshold)]["reslist"]



    def _export_structure(self, folder, filename, extension="cif", include_unused=True, include_misc=False) -> str|None:
        filename = "{}.structure.{}".format(filename, extension)
        filepath = os.path.join(folder, filename)
        self.paths["self"] = filepath
        if extension == "pdb":
            exp = bp.PDBIO()
            exp.set_structure(self)
            exp.save(filepath)
            return filepath
        elif extension == "cif":
            try:
                from .imports import write_atoms
                return write_atoms(self, filepath, include_unused=include_unused, include_misc=include_misc)
            except:
                return super()._export_structure(folder, filename, extension=extension, include_unused=include_unused, include_misc=include_misc)
        return None


    def show_exposed_residues(self, **kwargs):
        log("start", "swowing exposed residues")
        from src.bioiain.visualisation.pymol import PymolScript

        surfece_res_ids = self.get_surface_residues(**kwargs)
        print(surfece_res_ids)

        for atom in self.atoms():
            atom.set_bfactor(atom.get_misc("SASA"))
        path = self.export(folder="/tmp/bioiain/exports", structure=True, data=False, include_misc=True)

        threshold = kwargs.get("threshold", 50)
        script = PymolScript()
        script.load(path, self.get_name())
        script.color("all", "blue")
        for resn in self.data["surface"][str(threshold)]["reslist"]:
            script.color(f"i. {resn}", "orange")
        script.execute()
        log("end", "showing exposed residues")
