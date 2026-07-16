import sys

from src.bioiain.utilities import *

from src.bioiain.machine import FORCE

if "monomers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="monomers", oligo=0, ignore_blacklist=FORCE)

elif "multimers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.multimeric.list", name="multimers", oligo=1, ignore_blacklist=FORCE)

elif "pisa" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="pisa", oligo=0, ignore_blacklist=FORCE)
    dataset.add_list("./data/cath-dataset-nonredundant-S20.multimeric.list", oligo=1)

elif "receptors" in sys.argv:
    dataset = StructureDataset.from_list("./data/receptors.list", name="receptors", ignore_blacklist=FORCE)

elif "lbds" in sys.argv:
    dataset = StructureDataset.from_list("./data/lbds.list", name="lbds", ignore_blacklist=FORCE)

elif "aleph" in sys.argv:
    dataset = StructureDataset.from_list(["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"], name="aleph", ignore_blacklist=FORCE)

else:
    dataset = StructureDataset.from_list("./data/consensus.list", name="consensus", ignore_blacklist=FORCE)
