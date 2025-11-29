

import os
import sys

sys.path.append('..')

import src.bioiain as bi

bi.biopython.imports.read_mmcif("/home/iain/projects/vib-ai/internship/data/receptors/1M2Z.cif")