import sys
import Bio.PDB as bp

sys.path.append('.')
import src.bioiain as bi
from src.bioiain.utilities.logging import log
from src.bioiain.biopython import structure
from src.bioiain.symmetries.operations import *
from src.bioiain.symmetries.parsing import parse_crystal_card


pdb = parse_crystal_card("./test/1M2Z.pdb")
cif = parse_crystal_card("./test/1M2Z.cif")

print(pdb)
print(cif)





