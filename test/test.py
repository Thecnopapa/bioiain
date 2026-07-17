import os, json, sys, time, random
sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.base import *


log("start", "test.py")
log("title", "test.py")


from src.bioiain.tools.CCP4.TRACER import TRACER

example = {
    "cell": [10.51, 15.15, 6.54, 90, 151.7, 90],
    "d": 0.1,
    "centered": "C",
}
tracer = TRACER("test")
print(tracer)
tracer.run(**example)

exit()


entity = BIEntity.from_file("./5ezq.cif")
for res in entity.residues():
    if res.resnum == 477:
        print(res)
        for a in res.atoms:
            print(a)

entity.fragment()
exit()


term = logging.CursedTerminal()

heights = term.split_height(percentages=[10,50,40])
#print(heights)
widths = term.split_width(2)
#print(widths)

top =    term.add_window(heights[0], None,      heights[0], 0,         title=f"Top")
left =   term.add_window(heights[1], widths[0], heights[1], widths[0], title=f"Left")
rigth =  term.add_window(heights[1], widths[1], heights[1], widths[1], title=f"Right")
bottom = term.add_window(heights[2], None,      heights[2], 0,         title=f"Bottom")

top.print("data:")

top.print("aaaa")
top.print("bbb")
top.print("ccc")
top.print("ddd")
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

