import os, json, time, sys
from .logging import log





class SequenceMissmatchException(Exception):
	pass


class MisslabellingException(SequenceMissmatchException):
	pass

class DeletedIndex(Exception):
	def __init__(self, *args, next_=None, **kwargs):
		super().__init__(*args, **kwargs)
		self.next_n = next_n


class StructureLoadException(Exception):
	pass

class StructureRecoverException(StructureLoadException):
	pass