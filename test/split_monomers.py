import os, sys, json, shutil
sys.path.append('..')
from src.bioiain import log
from src.bioiain.utilities.exceptions import *
from src.bioiain.tools.PISA import *




log("start", "SPLITTING MONOMERS")
log("title", "SPLITTING MONOMERS")


from data import dataset


monf = f"./data/{dataset.name}.monomeric.list.temp"
multf = f"./data/{dataset.name}.multimeric.list.temp"
nopisaf = f"./data/{dataset.name}.nopisa.list.temp"


open(monf, "w")
open(multf, "w")
open(nopisaf, "w")


for entry in dataset:
    log("header", entry.name)
    code = entry.code
    fpath = entry.path
    file = os.path.basename(fpath)

    if os.path.getsize(fpath) > 1.5 * 1024 * 1024:
        log("warning", f"File: {file} too large! ({os.path.getsize(os.path.join(file_folder, file)) / 1024 / 1024 :3.2f} MiB)")
        with open(no_pisa_fname, "a") as f:
            f.write(code+"\n")
        continue

    pisa = PISA(pisa_id=code)
    try:
        pisa.analyse(fpath, force=False)
    except PISAError as e:
        log("warning", "PISA error at file:", file)
        log("warning", e)
        with open(nopisaf, "a") as f:
            f.write(code+"\n")
        continue

    if pisa["multimeric_state"] == 1:
        print(f" ... monomer found in {file}: {pisa['multimeric_state']}")
        with open(monf, "a") as f:
            f.write(code+"\n")

    else:
        print(f" ... multimer found in {file}: {pisa['multimeric_state']}")
        with open(multf, "a") as f:
            f.write(code+"\n")

    pisa.delete()

shutil.move(monf, monf.replace(".temp", ""))
shutil.move(multf, multf.replace(".temp", ""))
shutil.move(nopisaf, nopisaf.replace(".temp", ""))




log("title", "DONE")
log("end", "DONE")

