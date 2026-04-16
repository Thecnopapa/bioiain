
import os, json, sys


sys.path.append('..')
from src.bioiain.utilities import *
from src.bioiain.tools.PLINDER import PLINDERDatabase, PLINDERSystem


log("start", "plinder_test.py")
db = PLINDERDatabase(tutorial="--tutorial" in sys.argv)
print(db)

q = db.query(columns=["system_id", "entry_pdb_id", "entry_oligomeric_state"])
print(q)

system = PLINDERSystem(q.at[0,"system_id"])
print(system)
print(system.annotations())