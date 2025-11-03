import sys
import Bio.PDB as bp



sys.path.append('.')
import src.bioiain as bi
from src.bioiain.utilities.prints import log
from src.bioiain.biopython import structure
from src.bioiain.symmetries.operations import *

biop_structure = bp.Structure.Structure("4321")
test_structure = structure.Structure("1234")
print(test_structure)
print(type(test_structure))
print(test_structure, isinstance(test_structure,bp.Entity.Entity))
print(entity_to_frac(test_structure, {}))


print(biop_structure.__class__)
print(test_structure.__class__)

mix = structure.Structure.cast(biop_structure)
print(mix)
print(structure.Structure.__mro__)





