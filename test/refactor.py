import os, json, sys


sys.path.append('..')


from src.bioiain.base import *
from src.bioiain.base.mmcif import *
from src.bioiain.machine import *
from src.bioiain.machine.embeddings import CVEmbedding

from src.bioiain.visualisation.pymol import PymolScript


folder = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "1NOH", "6E52"], data_dir="./data" )

script = None
for file in os.listdir(folder):
    if "6E52" not in file:
        continue
    log("header", file)
    log("title", file)
    path = os.path.join(folder,  file)
    entity = FragmentedStructure.from_file(path).fragment()
    print(entity.ligands())

    #script = PymolScript()
    #script.load(path, entity.name())


    matrix = entity.cvmatrix()
    entity.export(target_folder = "trash")
    #print(matrix)
    #matrix.save_fig()
    entity.show_cvectors(execute=True)





