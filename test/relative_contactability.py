import os, sys, json
sys.path.append('..')

from src.bioiain.base import BIEntity
from src.bioiain.aleph import FragmentedStructure
from src.bioiain.machine.embeddings import ProteinEmbedding
from src.bioiain.machine.datasets import EmbeddingDataset
import vqvae_embeddings as embeddings
import polars as pl



def calculate_relative_contactability(dataset:EmbeddingDataset,
                                      contactability_pos:int=8,
                                      entity_class:type[BIEntity]=FragmentedStructure,
                                      embedding_class:type[ProteinEmbedding]=embeddings.ALEPHProteinEmbedding,
                                      padding:tuple[int, int]=(1,1)
                                      ) -> ProteinEmbedding:
    print(dataset)
    dataset.sequence_db(force=True)
    dataset.cluster(reassign=True, verbosity=3, force=True, linear=True)
    mmseqs = dataset.db()
    for e in dataset.embeddings.values():
        entity = entity_class.from_file(e["entity_path"], no_atoms=True)
        print(entity)
        query = entity.db().db_path()

        print(query, mmseqs)
        df = mmseqs.search(query)
        print(json.dumps(e, indent=4))
        embedding = embedding_class.from_file(e["embedding_path"])
        print(embedding)
        t = embedding.tensor()
        print(t)

        similar_ids = [i for i in df.get_column("target") if i != e["key"]]
        print("Similar IDs:", similar_ids)

        contactability = [float(r[contactability_pos]) for r in t]
        print(len(contactability))
        ps, pe = padding
        for key in similar_ids:
            row = df.row(by_predicate=pl.col("target") == key, named=True)
            print(row)
            target_embedding = embedding_class.from_file(dataset.embeddings[key]["embedding_path"])
            target_contactability = [float(r[contactability_pos]) for r in target_embedding.tensor()]
            print(len(target_contactability))
            qn = row["qstart"] - 1
            tn = row["tstart"] - 1
            for q, t in zip(row["qaln"][ps:-pe], row["taln"][ps:-pe]):
                print(qn, tn , q, t)
                addq, addt = True, True
                if q == "-":
                    addq=False
                if t == "-":
                    addt=False

                if addq and addt:
                    current_cont = contactability[qn]
                    extra_cont = target_contactability[tn]
                    print(qn, current_cont, tn, extra_cont)
                    contactability[qn] = contactability[qn] + extra_cont


                if addq:
                    qn += 1
                if addt:
                    tn += 1







        exit()


        # TODO: parse target ids
        # TODO: fetch embeddings for each target
        # TODO: parse aligned sequences
        # TODO: match sequence to embeddings
        # TODO: calculate relative contactability
        # TODO: generate new embeddings

        new_embedding = embedding_class.from_tensor(new_t)
        return new_embedding










if __name__ == "__main__":
    from src.bioiain.machine import EmbeddingDataset

    dataset = EmbeddingDataset(name=f"tokens_aleph_ExpandedALEPHEmbedding0")
    dataset.load()

    calculate_relative_contactability(dataset)

