import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .residue import Residue, DResidue
from ..utilities.sequences import d3to1
from ..utilities.logging import log


class Chain(bp.Chain.Chain, BiopythonOverlayClass):
    child_class = Residue
    disordered_child_class = DResidue

    def _init(self, *args, **kwargs):
        self.data["sequence"] = None


    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)


    def atoms(self, ca_only=False, hetatm=False, force=False, group_by_residue=False, disordered=False):
        from .imports import read_mmcif
        from .atom import BIAtom

        if self._atoms is None or force:
            #print("Reading atoms from CIF")
            atoms = read_mmcif(self.paths["self"], subset=["_atom_site"])("_atom_site")
            atoms = [BIAtom(a) for a in atoms]
            if not disordered:
                atoms = self._fix_disordered(atoms)
            self._atoms = atoms

        atoms = self._atoms
        if not hetatm:
            atoms = [a for a in atoms if a.type != "HETATM"]
        if ca_only:
            atoms = [a for a in atoms if a.name == "CA"]
        if group_by_residue:
            atoms_by_res = {}
            for atom in atoms:
                if atom.resnum in atoms_by_res:
                    atoms_by_res[atom.resnum].append(atom)
                else:
                    atoms_by_res[atom.resnum] = [atom]
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


    def claculate_sasa(self):
        pass


