import os, json, sys
sys.path.append('..')
from src.bioiain.aleph import FragmentedStructure


from src.bioiain.base import *
from src.bioiain.base.mmcif import *


folder = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf"], data_dir="./data" )

script = None
for file in os.listdir(folder):
    if "3BRF" not in file:
        continue
    log("header", file)
    log("title", file)
    path = os.path.join(folder,  file)
    entity = BIEntity.from_file(path)
    #entity = FragmentedStructure.from_file("./data/ccs/6F63.cif")
    entity = entity.fragment()


    entity.export()

    script = entity.show_cvectors(script = script, execute=False)

script.execute()