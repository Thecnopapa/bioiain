import os, json


from ..utilities.exceptions import *
from .atom import BIAtom



class BIResidue(object):
    child_class = BIAtom
    def __init__(self, atoms):
        if type(atoms) == dict:
            atoms = atoms.values()
        self.atoms = atoms
        self.ca = None
        self.cb = None
        self.c = None
        self.o = None
        self.n = None
        self.resnum = None
        self.resname = None
        self.resseq = None
        self.chain = None
        self.fragment = None
        self.is_residue = True

            

        for a in self.atoms:
            if len(a.resname) == 2:
                #log("Warning", "Nucleotide detected!, initialising BINucleotide instead")
                self.__class__ = BINucleoutide
                self.__init__(self.atoms)
                break

            if a.name == "CA":
                self.ca = a
            elif a.name == "CB":
                self.cb = a
            elif a.name == "C":
                self.c = a
            elif a.name == "O":
                self.o = a
            elif a.name == "N":
                self.n = a

        if self.is_residue:

            if len(self.atoms) == 1:
                log("Warning", "Only one atom given to residue, treating as CA")
                self.ca = self.atoms[0]
            if self.ca is None:
                log("error", "Trying to initialise residue with no CA")
                print([a.name for a in self.atoms])
                print()
                raise NoCaFound()

            self.fragment = self.ca.get_misc("fragment", None)

            self.resnum = self.ca.resnum
            self.resname = self.ca.resname
            self.resseq = self.ca.resseq
            self.chain = self.ca.chain
            if self.fragment is None:
                self.id = ( self.resname, self.resnum, self.resseq, self.chain)
            else:
                self.id = ( self.resname, self.resnum, self.resseq, self.chain, self.fragment)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} id={self.id}>"

    def to_atoms(self, key, value):
        for atom in self.atoms:
            atom.set_misc(key, value)

    def set_bfactor(self, bfactor):
        for a in self.atoms:
            a.set_bfactor(bfactor)
        return self

class BINucleoutide(object):
    child_class = BIAtom
    def __init__(self, atoms):
        if type(atoms) == dict:
            atoms = atoms.values()
        self.atoms = atoms
        self.c1 = None
        self.resnum = None
        self.resname = None
        self.resseq = None
        self.chain = None
        self.fragment = None
        self.is_residue = False

        for a in self.atoms:
            #print(a)
            #print(a.name)
            if a.name == "C1":
                self.c1 = a
                break

        if self.c1 is None:
            log("error", "Trying to initialise nucleotide with no C1")
            print([a.name for a in self.atoms])
            print()
            raise NoCaFound()

        self.fragment = self.c1.get_misc("fragment", None)

        self.resnum = self.c1.resnum
        self.resname = self.c1.resname
        self.resseq = self.c1.resseq
        self.chain = self.c1.chain
        if self.fragment is None:
            self.id = ( self.resname, self.resnum, self.resseq, self.chain)
        else:
            self.id = ( self.resname, self.resnum, self.resseq, self.chain, self.fragment)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__} id={self.id}>"

    def to_atoms(self, key, value):
        for atom in self.atoms:
            atom.set_misc(key, value)
