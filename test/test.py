import os, json, sys

import torch

sys.path.append('..')

from src.bioiain.utilities import *

log("start", "test.py")
log("title", "test.py")



from src.bioiain.utilities.logging import *
from src.bioiain.base import *
#from compactness_base import *
#from compactness_models import *
#from src.bioiain.machine import *
#set_seed()



entity = BIEntity.from_file("2plj.cif", export_folder="trash")
print(entity, entity.path())

entity.img3D(property="b")


