import os, json
import numpy as np


import Bio.PDB as bp
from .base import BiopythonOverlayClass
from .atom import Atom, DAtom, BIAtom

from ..utilities.logging import log
from ..utilities.exceptions import *



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
            print(a)
            print(a.name)
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






class Residue(bp.Residue.Residue, BiopythonOverlayClass):
    child_class = Atom
    disordered_child_class = DAtom

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

class DResidue(bp.Residue.DisorderedResidue, Residue):
    child_class = Atom
    disordered_child_class = DAtom

    def __repr__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __str__(self):
        return "<bi.{} id={}>".format(self.__class__.__name__, self.id)

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        for a in self.disordered_get_id_list():
            self[a] = Residue.cast(self.disordered_get(a))
        self.disordered_select(self.id)