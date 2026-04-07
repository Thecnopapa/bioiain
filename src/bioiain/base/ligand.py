import os, json

from ..utilities.exceptions import *
from ..utilities import *
from .atom import BIAtom








class Water(object):
	child_class = BIAtom
	type = "water"
	def __init__(self, atoms, **kwargs):
		self.atoms = atoms

		for a in self.atoms:
			if a.name == "O":
				self.o = a
		self.resseq = self.o.resseq
		self.id = self.resseq
		self.relevant = False

	def __repr__(self):
		return f"<bi.{self.__class__.__name__} id={self.id}>"


class Ligand(object):
	child_class = BIAtom
	type = "ligand"
	def __init__(self, atoms, parent=None, **kwargs):
		self.atoms = atoms

		self.name = self.atoms[0].resname
		self.chain = self.atoms[0].chain
		self.complex = self.atoms[0].complex
		self.entity = self.atoms[0].entity
		self.model = self.atoms[0].model

		self.id = (self.name, self.complex)
		self.id2 = (self.name, self.chain)

		self.relevant = False

		self._com = None

		if parent is not None:
			self._determine_relevance(entity=parent)


	def _determine_relevance(self, entity=None):
		#if self.name == "DEX":
		#	self.relevant = True
		if entity is None:
			self.relevant = True
		else:
			sa = self._calculate_sasa(entity=entity)
			print("SA", sa)

		return self.relevant

	def _calculate_sasa(self, entity=None):
		from ..tools.SASA import SASA
		sasa = SASA()
		print(sasa.compute(entity=entity, targets=self.atoms))
		print(sasa)
		exit()



	def com(self, force=False):
		if self._com is None or force:
			from ..utilities.maths import find_com
			self._com = find_com(self.atoms)
		return self._com

	def __repr__(self):
		return f"<bi.{self.__class__.__name__} id={self.id2}({self.complex})>"