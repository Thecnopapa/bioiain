import os, sys, json
sys.path.append('..')



def calculate_relative_contactability(dataset):
    print(dataset)
    dataset.sequence_db(force=True)
    dataset.cluster(reassign=True, verbosity=3, force=True)








if __name__ == "__main__":
    from src.bioiain.machine import EmbeddingDataset

    dataset = EmbeddingDataset(name=f"tokens_aleph_ExpandedALEPHEmbedding0")
    dataset.load()

    calculate_relative_contactability(dataset)

