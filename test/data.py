import sys

from src.bioiain.utilities import *



if "monomers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="monomers", oligo=0)

elif "multimers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.multimeric.list", name="multimers", oligo=1)

elif "pisa" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="pisa", oligo=0)
    dataset.add_list("./data/cath-dataset-nonredundant-S20.multimeric.list", oligo=1)

elif "receptors" in sys.argv:
    dataset = StructureDataset.from_list("./data/receptors.list", name="receptors")

elif "lbds" in sys.argv:
    dataset = StructureDataset.from_list("./data/lbds.list", name="lbds")

elif "aleph" in sys.argv:
    dataset = StructureDataset.from_list(["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"], name="aleph")

else:
    dataset = StructureDataset.from_list("./data/consensus.list", name="consensus")
