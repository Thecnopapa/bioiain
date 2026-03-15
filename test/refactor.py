import os, json, sys
sys.path.append('..')
from src.bioiain.aleph import FragmentedStructure


from src.bioiain.base import *
from src.bioiain.base.mmcif import *

entity = FragmentedStructure.from_file("./3HHB.cif")

entity = entity.fragment_with_aleph()




entity.export()

entity.show_fragments()