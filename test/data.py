import sys

from src.bioiain.utilities import *


ignore_blacklist = ("--force" in sys.argv) or ("-f" in sys.argv) or ("--no-blacklist" in sys.argv)


if "monomers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="monomers", oligo=0, ignore_blacklist=ignore_blacklist)

elif "multimers" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.multimeric.list", name="multimers", oligo=1, ignore_blacklist=ignore_blacklist)

elif "pisa" in sys.argv:
    dataset = StructureDataset.from_list("./data/cath-dataset-nonredundant-S20.monomeric.list", name="pisa", oligo=0, ignore_blacklist=ignore_blacklist)
    dataset.add_list("./data/cath-dataset-nonredundant-S20.multimeric.list", oligo=1)

elif "receptors" in sys.argv:
    dataset = StructureDataset.from_list("./data/receptors.list", name="receptors", ignore_blacklist=ignore_blacklist)

elif "lbds" in sys.argv:
    dataset = StructureDataset.from_list("./data/lbds.list", name="lbds", ignore_blacklist=ignore_blacklist)

elif "aleph" in sys.argv:
    dataset = StructureDataset.from_list(["1M2Z", "3HBB", "6F63", "5LXN", "3brf", "6e52", "7t2y", "3kg2", "2GEJ", "2bis"], name="aleph", ignore_blacklist=ignore_blacklist)

elif "ccs" in sys.argv:
        dataset = StructureDataset.from_list("./data/ccs.list", name="ccs", ignore_blacklist=ignore_blacklist)

elif "test" in sys.argv:
        dataset = StructureDataset.from_list(["1A92"], name="test", ignore_blacklist=ignore_blacklist)

else:
    dataset = StructureDataset.from_list("./data/consensus.list", name="consensus", ignore_blacklist=ignore_blacklist)
