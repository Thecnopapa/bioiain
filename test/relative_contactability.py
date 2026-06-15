import os, sys, json


sys.path.append('..')

from src.bioiain.base import BIEntity
from src.bioiain.aleph import FragmentedStructure
from src.bioiain.machine.embeddings import ProteinEmbedding
from src.bioiain.machine.datasets import EmbeddingDataset
import vqvae_embeddings as embeddings
import polars as pl
import numpy as np
from torch import Tensor



def calculate_relative_contactability(dataset:EmbeddingDataset,
                                      contactability_pos:int=8,
                                      entity_class:type[BIEntity]=FragmentedStructure,
                                      embedding_class:type[ProteinEmbedding]=embeddings.ALEPHProteinEmbedding,
                                      padding:tuple[int, int]=(1,1)
                                      ) -> EmbeddingDataset:
    print(dataset)
    print(dataset.data["embedding_class"], embedding_class.__name__)
    assert dataset.data["embedding_class"] == embedding_class.__name__
    dataset.sequence_db(force=True)
    dataset.cluster(reassign=True, verbosity=3, force=True, linear=True)
    mmseqs = dataset.db()
    new_dataset = dataset.__class__(name=dataset.data["name"]+"_relative")
    for e in dataset.embeddings.values():
        entity = entity_class.from_file(e["entity_path"], no_atoms=True)
        print(entity)
        print(len(entity.cvectors()))
        query = entity.db().db_path()

        print(query, mmseqs)
        df = mmseqs.search(query)
        print(json.dumps(e, indent=4))
        embedding = embedding_class.from_file(e["embedding_path"])
        print(embedding)
        t = embedding.tensor()
        print(t)
        print(t.shape)
        print("####")

        similar_ids = [i for i in df.get_column("target") if i != e["key"]]
        print("Similar IDs:", similar_ids)
        print(t)
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

                #print(qn, tn , q, t)
                addq, addt = True, True
                if q == "-":
                    addq=False
                if t == "-":
                    addt=False

                if addq and addt:
                    extra_cont = target_contactability[tn]
                    #print(qn, current_cont, tn, extra_cont)
                    contactability[qn] = contactability[qn] + extra_cont


                if addq:
                    qn += 1
                if addt:
                    tn += 1


            print(len(contactability))
            print(contactability)

        if len(similar_ids) > 0:
            contactability = np.array(contactability)
            contactability = contactability / (len(similar_ids)+1)
            print(len(contactability))
            print(contactability)

            new_tensor = np.array(t)
            for r, c in zip(new_tensor, contactability):
                print(r, c)
                r[contactability_pos] = c
            new_tensor = Tensor(new_tensor)
            new_embedding = embeddings.RelativeALEPHEmbedding.from_tensor(new_tensor)


        else:
            new_embedding = embeddings.RelativeALEPHEmbedding.from_tensor(t)
        new_embedding.entity_path = embedding.entity_path
        new_embedding.sequence = embedding.sequence
        new_embedding.entity = embedding.entity
        new_embedding.param_names = embedding.param_names
        new_embedding.missing_indexes = embedding.missing_indexes
        new_embedding.residue_embedding_class = embedding.residue_embedding_class


        new_embedding.save()
        print(new_embedding)
        new_dataset.add(new_embedding)
    return new_dataset











if __name__ == "__main__":
    from src.bioiain.machine import EmbeddingDataset

    dataset = EmbeddingDataset(name=f"tokens_aleph_ExpandedALEPHEmbedding0")
    dataset.load()

    calculate_relative_contactability(dataset)

