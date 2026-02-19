import os, json, time, sys
from .logging import log





class SequenceMissmatchException(Exception):
	pass


class MisslabellingException(SequenceMissmatchException):
	pass

class DeletedIndex(Exception):
	pass


class StructureLoadException(Exception):
	pass

class StructureRecoverException(StructureLoadException):
	pass