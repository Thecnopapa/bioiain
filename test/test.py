import os, json, sys, time
sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.base import *


log("start", "test.py")
log("title", "test.py")







progress = logging.GraphProgress()

import random
for n in range(50):
	progress.update(random.randrange(1,20))
	time.sleep(1)

progress.end()
exit()




from data import dataset

from saprot3D import *
from compactness_base import *
from compactness_models import *

IMG_SIZE = 16



for entity in dataset.entities(entity_class=CompactStructure):

	entity.img3D(property="compactness", show_plot=False, gif=True)

