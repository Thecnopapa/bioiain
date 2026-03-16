import os, json, sys
sys.path.append('..')
from src.bioiain.aleph import FragmentedStructure


from src.bioiain.base import *
from src.bioiain.base.mmcif import *

entity = BIEntity.from_file("./3HHB.cif")
#entity = FragmentedStructure.from_file("./data/ccs/6F63.cif")
entity = entity.fragment()




entity.export()

entity.show_cvectors()