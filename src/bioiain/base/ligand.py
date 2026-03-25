import os, json

from ..utilities.exceptions import *
from ..utilities import *
from .atom import BIAtom








class Water(object):
	child_class = BIAtom
	type = "water"
	def __init__(self, atoms):
		self.atoms = atoms

		for a in self.atoms:
			if a.name == "O":
				self.o = a
		self.resseq = self.o.resseq
		self.id = self.resseq

	def __repr__(self):
		return f"<bi.{self.__class__.__name__} id={self.id}>"


class Ligand(object):
	child_class = BIAtom
	type = "ligand"
	def __init__(self, atoms):
		self.atoms = atoms

		self.name = self.atoms[0].resname
		self.chain = self.atoms[0].chain
		self.complex = self.atoms[0].complex
		self.entity = self.atoms[0].entity
		self.model = self.atoms[0].model

		self.id = (self.name, self.complex)
		self.id2 = (self.name, self.chain)




	def __repr__(self):
		return f"<bi.{self.__class__.__name__} id={self.id2}({self.complex})>"