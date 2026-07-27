import os, json, sys
sys.path.append('..')


from src.bioiain.utilities import MMSEQS2



msa = MMSEQS2("/localdata/iain/bioiain/test/bioiain.d/datasets/tokens_aleph_v1.dataset.fasta")
msa.cluster()




