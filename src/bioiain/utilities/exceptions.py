import os, json, time, sys
from .logging import log




# Labelling related
class SequenceMissmatchException(Exception):
	pass


class MisslabellingException(SequenceMissmatchException):
	pass


# Dataset related
class DeletedIndex(Exception):
	def __init__(self, *args, next_n=None, **kwargs):
		super().__init__(*args, **kwargs)
		self.next_n = next_n


# MMCIF related
class MMCIFError(Exception):
	pass
class MMCIFTypeError(MMCIFError):
	pass


# Structure import related
class StructureLoadException(Exception):
	pass

class StructureRecoverException(StructureLoadException):
	pass
class StructureNotFound(StructureRecoverException):
	pass

class AlreadyLoaded(StructureLoadException):
	pass
class UnknownFormat(StructureLoadException):
	pass


# CCP4 related
class CCP4Error(Exception):
	pass

class CCP4NotEnabled(CCP4Error):
	pass

class PISAError(CCP4Error):
	pass


# Model related
class ModelNotFound(Exception):
	pass

# Residue related
class NoCaFound(Exception):
	pass

class NoBackbone(Exception):
	pass
	
class NoMatchingClass(Exception):
	pass



# Crystal related
class FractionalConversionError(Exception):
	pass

class AlreadyFractional(FractionalConversionError):
	pass

class AlreadyOrthogonal(FractionalConversionError):
	pass

class CrystalError(Exception):
	pass

class MissingCrystalInfo(CrystalError):
	pass

# ALEPH related

class ALEPHError(Exception):
	pass


# CVector Related
class CVMatrixError(Exception):
	pass
class NoNeighboursFound(CVMatrixError):
	pass

# PLINDER related

class PLINDERError(Exception):
	pass

class PLINDERSystemNotLoaded(PLINDERError):
	pass


# MMSEQS2 related

class MMSEQS2Error(Exception):
	pass

class DatabaseError(MMSEQS2Error):
	pass

class ClusteringError(MMSEQS2Error):
	pass


# Embedding related
class NoEmbeddingForThisResidue(Exception):
	pass



















