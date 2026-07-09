import os, json, sys, time, random
sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.base import *


log("start", "test.py")
log("title", "test.py")







term, top, left, right, bottom = logging.quad_term()
top.print("data:")

top.print("aaaa")
top.print("bbb", end = "\r")
top.print("ccc")
top.print("ddd", end=" ")
top.print("eee")
left.print("default", c="default")
left.print("black", c="black")
left.print("blue", c="blue")
left.print("cyan", c="cyan")
left.print("green", c="green")
left.print("magenta", c="magenta")
left.print("red", c="red")
left.print("white", c="white")
left.print("yellow", c="yellow")


time.sleep(2)
#print(term)

time.sleep(10)
term.close()




exit()




from data import dataset

from saprot3D import *
from compactness_base import *
from compactness_models import *

IMG_SIZE = 16



for entity in dataset.entities(entity_class=CompactStructure):

    entity.img3D(property="compactness", show_plot=False, gif=True)

