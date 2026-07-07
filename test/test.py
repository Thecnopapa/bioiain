import os, json, sys
sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.base import *


log("start", "test.py")
log("title", "test.py")

from data import dataset



from saprot3D import *
from compactness_base import *
from compactness_models import *

IMG_SIZE = 16


entity = next(dataset.entities())

entity.img3D(property="b", plot=True)

exit()
