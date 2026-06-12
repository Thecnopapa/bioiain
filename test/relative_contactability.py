import os, sys, json
sys.path.append('..')

from src.bioiain.aleph import FragmentedStructure




def calculate_relative_contactability(dataset):
    print(dataset)
    dataset.sequence_db(force=True)
    dataset.cluster(reassign=True, verbosity=3, force=True, linear=True)
    mmseqs = dataset.db()
    for e in dataset.embeddings.values():
        entity = FragmentedStructure.from_file(e["entity_path"], no_atoms=True)
        print(entity)
        query = entity.db().db_path()

        print(query, mmseqs)
        mmseqs.search(query)









if __name__ == "__main__":
    from src.bioiain.machine import EmbeddingDataset

    dataset = EmbeddingDataset(name=f"tokens_aleph_ExpandedALEPHEmbedding0")
    dataset.load()

    calculate_relative_contactability(dataset)

