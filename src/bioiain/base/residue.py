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

        for a in self.atoms:
            #print(a)
            #print(a.name)
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

        if len(self.atoms) == 1:
            log("Warning", "Only one atom given to residue, treating as CA")
            self.ca = self.atoms[0]
        if self.ca is None:
            log("error", "Trying to initialise residue with no CA")
            raise NoCaFound()

        self.resnum = self.ca.resnum
        self.resname = self.ca.resname
        self.resseq = self.ca.resseq
        self.chain = self.ca.chain
        self.id = ( self.resname, self.resnum, self.resseq, self.chain)

    def __repr__(self):
        return f"<bi.BIResidue id={self.id}>"

