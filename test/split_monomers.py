import os, sys, json
sys.path.append('..')
from src.bioiain.biopython import downloadPDB
from src.bioiain import log
from src.bioiain.utilities.parallel import *
from src.bioiain.utilities.exceptions import *
from src.bioiain.tools.PISA import *




log("start", "SPLITTING MONOMERS")
log("title", "SPLITTING MONOMERS")


file_folder = downloadPDB("./data", "cath-nonredundant-S20",
                                               file_path="./data/cath-dataset-nonredundant-S20.list",
                                               file_format="cif",
                                               overwrite=False) 

monomeric_fname = "./data/cath-dataset-nonredundant-S20.monomeric.list"
multimeric_fname = "./data/cath-dataset-nonredundant-S20.multimeric.list"
no_pisa_fname = "./data/cath-dataset-nonredundant-S20.nopisa.list"


monf = open(monomeric_fname, "w")
multf = open(multimeric_fname, "w")
nopisaf = open(no_pisa_fname, "w")


for file in sorted(os.listdir(file_folder)):
	log("header", file)
	code = file.split(".")[0]
	fpath = os.path.join(file_folder, file)

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





















log("title", "DONE")
log("end", "DONE")

