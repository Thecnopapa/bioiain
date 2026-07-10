import os, sys, json
sys.path.append('..')
from src.bioiain import log
from src.bioiain.utilities.exceptions import *
from src.bioiain.tools.PISA import *




log("start", "SPLITTING MONOMERS")
log("title", "SPLITTING MONOMERS")


from data import dataset


monf = open(f"{dataset.name}.monomeric.list.temp", "w")
multf = open(f"{dataset.name}.multimeric.list.temp", "w")
nopisaf = open(f"{dataset.name}.nopisa.list.temp", "w")


for entry in dataset:
    log("header", entry["name"])
    code = entry["code"]
    fpath = entry["path"]

    if os.path.getsize(os.path.join(file_folder, file)) > 1.5 * 1024 * 1024:
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
        with open(no_pisa_fname, "a") as f:
            f.write(code+"\n")
        continue

    if pisa["multimeric_state"] == 1:
        print(f" ... monomer found in {file}: {pisa['multimeric_state']}")
        with open(monomeric_fname, "a") as f:
            f.write(code+"\n")

    else:
        print(f" ... multimer found in {file}: {pisa['multimeric_state']}")
        with open(multimeric_fname, "a") as f:
            f.write(code+"\n")

    pisa.delete()

os.copy(monomeric_fname, monomeric_fname.replace(".temp", ""))
os.copy(multimeric_fname, multimeric_fname.replace(".temp", ""))
os.copy(no_pisa_fname, no_pisa_fname.replace(".temp", ""))

os.remove(monomeric_fname)
os.remove(multimeric_fname)
os.remove(no_pisa_fname)


log("title", "DONE")
log("end", "DONE")

