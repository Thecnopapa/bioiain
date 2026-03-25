import os, json, sys


sys.path.append('..')


from src.bioiain.base import *
from src.bioiain.base.mmcif import *
from src.bioiain.machine import *
from src.bioiain.machine.embeddings import CVEmbedding

from src.bioiain.visualisation.pymol import PymolScript


folder = downloadPDB( list_name="aleph", pdb_list=["1M2Z", "3HBB", "6F63", "5LXN", "3brf"], data_dir="./data" )

script = None
for file in os.listdir(folder):
    if "1M2Z" not in file:
        continue
    log("header", file)
    log("title", file)
    path = os.path.join(folder,  file)
    entity = BIEntity.from_file(path)

    script = PymolScript()
    script.load(path, entity.name())

    for a in entity.atoms():
        #print(entity.symops())
        if a.name != "CA":
            continue
        for opn, symop in entity.symops(as_dict=True).items():
            #print("OPN", opn)
            if opn != 2:
                continue
            symcoords = a.at(symop, entity.params())
            for symcoord in symcoords:
                script.line(name=f"op{opn}", coord1=a.coord, coord2=symcoord)

    for opn in entity.symops():
        script.group(f"op{opn}")

    script.execute()


