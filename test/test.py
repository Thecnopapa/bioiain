import os, json, sys

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")

from src.bioiain.utilities.logging import *
from src.bioiain.base import *
from compactness import *

entity = CompactStructure.from_file(os.path.join(".", "1M2Z.cif"), code="TEST", force=False)
print(entity)


entity._calculate_compactness()







exit()




