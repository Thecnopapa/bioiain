import sys
import Bio.PDB as bp

sys.path.append('.')
import src.bioiain as bi
from src.bioiain.utilities.logging import log
from src.bioiain.biopython import imports
from src.bioiain.symmetries.operations import *
from src.bioiain.symmetries.parsing import parse_crystal_card
from src.bioiain.utilities.entities import print_all_coords




pdb = imports.loadPDB("./test/1M2Z.pdb")
cif = imports.loadPDB("./test/1M2Z.cif")




print(pdb.__dict__.keys())
print(cif.__dict__.keys())




