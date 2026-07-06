
from ..utilities.exceptions import *
from .atom import BIAtom
from .ligand import Ligand, Water
from ..utilities import d3to1





def build_res(atoms, ignore_errors=False, **kwargs):
    try:
        resnames = [a.resname for a in atoms]
        atomnames = [a.name for a in atoms]

        if resnames == ["HOH"]:
            return Water(atoms, **kwargs)

        if max(len(r) for r in resnames) < 3:
            return BINucleotide(atoms, **kwargs)


        if "CA" in atomnames:
            return BIResidue(atoms, **kwargs)

        elif "C6" in atomnames:
            return BIHexose(atoms, **kwargs)
        elif "C5" in atomnames:
            return BIRibose(atoms, **kwargs)

        atomtypes = [a.type for a in atoms]
        if atomtypes == ["HETATM"]:
            return Ligand(atoms, **kwargs)

        raise NoMatchingClass()

    except (NoMainAtomFound, NoBackbone, NoMatchingClass, NotImplementedError) as e:
        if ignore_errors:
            log("warning", e)
            return None
        else:
            raise







class BIResidue(object):
    child_class = BIAtom
    type="residue"
    main_atom_name = "CA"
    def __init__(self, atoms, require_ca=True, **kwargs):
        if type(atoms) == dict:
            atoms = atoms.values()
        self.atoms = [a for a in atoms if a.element != "H"]
        self.ca = None
        self.cb = None
        self.c = None
        self.o = None
        self.n = None
        self.resnum = None
        self.resname = None
        self.rn1 = None
        self.resseq = None
        self.chain = None
        self.complex = None
        self.fragment = None
        self.is_residue = True
        self.is_disordered = False

        self._sasa = None
        self._av_sasa = None
        self._norm_sasa = None


            

        for a in self.atoms:

            if a.name == self.main_atom_name:
                self.ca = a
            elif a.name == "CB":
                self.cb = a
            elif a.name == "C":
                self.c = a
            elif a.name == "O":
                self.o = a
            elif a.name == "N":
                self.n = a
        self.backbone = [self.ca, self.c, self.o, self.n]

        if self.is_residue:

            if len(self.atoms) == 1:
                log("Warning", "Only one atom given to residue, treating as CA")
                self.ca = self.atoms[0]
            if self.ca is None:
                if require_ca:
                    print("No CA")
                    print([a.name for a in self.atoms])
                    print([a.atomnum for a in self.atoms])
                    raise NoCaFound("Trying to initialise residue with no CA")

            self.set_fragment()

            self.resnum = self.ca.resnum
            self.resname = self.ca.resname
            try:
                self.rn1 = d3to1[self.resname]
            except:
                self.rn1 = "X"

            self.resseq = self.ca.resseq
            self.chain = self.ca.chain
            self.entity = self.ca.entity
            self.complex = self.ca.complex

            self.is_disordered = not self.ca.ins_code is None
            if self.is_disordered:
                raise NotImplementedError()


            if any([a is None for a in self.backbone]):
                print("No backbone")
                print(self)
                print(self.backbone)
                raise NoBackbone(self, f"Trying to initialise {self.__class__.__name__} with no backbone")

    def id(self): return self.resname, self.resnum, self.resseq, self.chain, self.entity, self.complex, self.fragment

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} id={self.id()}>"

    def name(self):
        return "_".join([str(v) for v in self.id()])

    def to_atoms(self, key, value):
        for atom in self.atoms:
            atom.set_misc(key, value)

    def bfactor(self):
        return self.ca.b

    def set_misc(self, key, value):
        for a in self.atoms:
            a.set_misc(key, value)
        return self

    def get_misc(self, key, keyerror="undefined"):
        return self.ca.get_misc(key, keyerror=keyerror)

    def set_bfactor(self, bfactor):
        for a in self.atoms:
            a.set_bfactor(bfactor)
        return self

    def set_fragment(self, fragment=None, unset=False):
        if fragment is not None or unset:
            for a in self.atoms:
                a.set_misc(fragment)
        self.fragment = self.ca.get_misc("fragment", None)

    def sasa(self, normalised=False, average=False, force=False):
        if normalised:
            if self._norm_sasa is None or force:
                self._read_sasa()
            return self._norm_sasa
        elif average:
            if self._av_sasa is None or force:
                self._read_sasa()
            return self._av_sasa
        else:
            if self._sasa is None or force:
                self._read_sasa()
            return self._sasa

    def _read_sasa(self):
        from ..tools.SASA import residue_sasas
        summ = 0
        total = 0
        for a in self.atoms:
            s = a.get_misc("SASA")
            summ += s
            total += 1
        av = summ / total if total != 0 else None
        self._sasa = summ
        self._av_sasa = av
        self._norm_sasa = min(1, summ / residue_sasas[self.resname]) if self.resname in residue_sasas and av is not None else None
        #print( self._av_sasa, self._sasa, residue_sasas[self.resname], self._norm_sasa)
        return self._sasa





class BISugar(object):
    child_class = BIAtom
    type="ligand"
    main_atom_name="C1"
    def __init__(self, atoms, **kwargs):
        if type(atoms) == dict:
            atoms = atoms.values()
        self.atoms = atoms
        self.main = None
        self.resnum = None
        self.resname = None
        self.resseq = None
        self.chain = None
        self.fragment = None
        self.is_residue = False



        for a in self.atoms:
            if a.name == self.main_atom_name:
                self.main = a
                break

        if self.main is None:
            print("No main", self.__class__.__name__)
            print([a.name for a in self.atoms])
            print([a.atomnum for a in self.atoms])
            [print(a, a.id4()) for a in self.atoms]
            raise NoMainAtomFound(f"Trying to initialise {self.__class__.__name__} with no {self.main_atom_name}")

        self.fragment = self.main.get_misc("fragment", None)

        self.resnum = self.main.resnum
        self.resname = self.main.resname
        self.resseq = self.main.resseq
        self.chain = self.main.chain


    def id(self): return self.resname, self.resnum, self.resseq, self.chain, self.fragment

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} id={self.id()}>"

    def to_atoms(self, key, value):
        for atom in self.atoms:
            atom.set_misc(key, value)


class BIRibose(BISugar):
    pass

class BIHexose(BISugar):
    pass

class BINucleotide(BIHexose):
    type="nucleotide"
    main_atom_name="N1"

class BIRNA(BINucleotide):
    pass

class BIDBNA(BINucleotide):
    pass
