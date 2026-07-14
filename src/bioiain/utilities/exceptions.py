import os, json, time, sys
from .logging import log


# Miscellaneous
class NotAGoodIdea(Exception):
    pass


# 3rd party related
class MissingProgram(Exception):
    pass



#PDB download related
class DownloadError(Exception):
    pass


# Labelling related
class SequenceMissmatchException(Exception):
    pass

class MisslabellingException(SequenceMissmatchException):
    pass

class ChainMissmatchException(Exception):
    pass

class MultipleChainsDetected(ChainMissmatchException):
    pass

class NoChainsDetected(ChainMissmatchException):
    pass


# EmbeddingDataset related
class DeletedIndex(Exception):
    def __init__(self, *args, next_n=None, **kwargs):
        super().__init__(*args, **kwargs)
        self.next_n = next_n

class EmbeddingDatasetNotFound(FileNotFoundError):
    pass

class AlreadyInDataset(Exception):
    pass



# MMCIF related
class MMCIFError(Exception):
    pass
class MMCIFTypeError(MMCIFError):
    pass


# Sequence related
class SequenceNotFound(Exception):
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

class DisorderParsingError(StructureLoadException):
    pass

# Structure export related
class StructureExportError(Exception):
    pass
class ExportingNoAtoms(StructureExportError):
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
class ModelSaveError(Exception):
    pass
class TryingToSaveInferenceModel(ModelSaveError):
    pass


# LossRelated
class LossIsZero(Exception):
    pass


# Residue related

class ResidueBuildingError(Exception):
    def __init__(self, residue, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.residue = residue
    
class NoMainAtomFound(ResidueBuildingError):
    pass

class NoCaFound(NoMainAtomFound):
    pass

class NoBackbone(ResidueBuildingError):
    pass

class NoMatchingClass(Exception):
    pass

# Nucleotide related



# Crystal related
class FractionalConversionError(Exception):
    pass

class AlreadyFractional(FractionalConversionError):
    pass

class AlreadyOrthogonal(FractionalConversionError):
    pass

class CrystalError(Exception):
    pass

class MissingCrystalInfo(CrystalError, StructureLoadException):
    pass

# ALEPH related

class ALEPHError(Exception):
    pass

class ALEPHMissmatch(ALEPHError):
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

class SearchError(MMSEQS2Error):
    pass

class TsvError(MMSEQS2Error):
    pass

# Foldseek related

class FoldseekError(Exception):
    pass

class FoldseekParsingError(FoldseekError):
    pass

class NameParsingError(FoldseekParsingError):
    pass


# Embedding related
class EmbeddingGenerationError(Exception):
    pass

class NoEmbeddingForThisResidue(EmbeddingGenerationError):
    pass

class NoEmbeddingForThisProtein(EmbeddingGenerationError):
    pass

class EmptyTensor(Exception):
    pass

class EmbeddingLoadError(Exception):
    pass

class NoTensorAvailable(EmbeddingLoadError):
    pass

class UnknownEmbeddingFormat(EmbeddingLoadError):
    pass
