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
from src.bioiain.utilities.exceptions import *
from src.bioiain import log










class ContactactStructure(BIEntity):
    pass
























class RelativeContactabilityError(Exception):
    pass


def calculate_relative_contactability(dataset:EmbeddingDataset,
                                      contactability_pos:int=8,
                                      entity_class:type[BIEntity]=FragmentedStructure,
                                      embedding_class:type[ProteinEmbedding]=embeddings.ALEPHProteinEmbedding,
                                      ) -> EmbeddingDataset:
    print(dataset)
    print(dataset.data["embedding_class"], embedding_class.__name__)
    assert dataset.data["embedding_class"] == embedding_class.__name__
    dataset.sequence_db(force=False)
    dataset.cluster(reassign=True, force=False, linear=True)
    mmseqs = dataset.db()
    new_dataset = dataset.__class__(name=dataset.data["name"]+"_relative")
    for e in dataset.embeddings.values():
        try:
            entity = entity_class.from_file(e["entity_path"], no_atoms=True)
            #print(entity)
            #print(len(entity.cvectors()))
            try:
                query = entity.db().db_path()
            except SequenceNotFound:
                entity = entity_class.from_file(e["entity_path"], no_atoms=False)
                try:
                    query = entity.db().db_path()
                except:
                    log("warning", "No sequence for:", entity)
                    continue



            #print(query, mmseqs)
            df = mmseqs.search(query, dataset_name=str(dataset))
            #print(json.dumps(e, indent=4))
            embedding = embedding_class.from_file(e["embedding_data"])
            #print(embedding)
            tensor = embedding.tensor()
            #print(tensor)
            print(tensor.shape)
            #print("####")

            similar_ids = [i for i in df.get_column("target") if i != e["key"]]
            print("Similar IDs:", similar_ids)
            #print(tensor)
            contactability = [float(r[contactability_pos]) for r in tensor]
            #print(len(contactability))
            for key in similar_ids:
                row = df.row(by_predicate=pl.col("target") == key, named=True)
                #print(row)
                target_embedding = embedding_class.from_file(dataset.embeddings[key]["embedding_data"])
                target_contactability = [float(r[contactability_pos]) for r in target_embedding.tensor()]
                #print(len(target_contactability))
                #qn, tn = 0, 0
                qn = row["qstart"] - 1
                tn = row["tstart"] - 1
                #print(len(row["qaln"]), len(row["taln"]))
                for q, t in zip(row["qaln"], row["taln"]):
                    #print(qn, tn , q, t)

                    addq, addt = True, True
                    skip = False

                    if qn in embedding.missing_indexes:
                        skip = True
                        addq = False
                    if tn in target_embedding.missing_indexes:
                        skip = True
                        addt = False

                    if q == "-":
                        addq=False
                    if t == "-":
                        addt=False

                    if addq and addt and (not skip):
                        extra_cont = target_contactability[tn]
                        print(qn, len(contactability), tn, len(target_contactability), end="\r")
                        try:
                            contactability[qn] = contactability[qn] + extra_cont
                        except IndexError:
                            print()
                            print("qn/tn", qn, tn)
                            print("seq_len", len(embedding.sequence), len(target_embedding.sequence))
                            print("seq", embedding.sequence)
                            print("missing", embedding.missing_indexes, target_embedding.missing_indexes)
                            print(qn, len(contactability))
                            raise RelativeContactabilityError()
                    if addq:
                        qn += 1
                    if addt:
                        tn += 1


            if len(similar_ids) > 0:
                contactability = np.array(contactability)
                contactability = contactability / (len(similar_ids)+1)
                #print(len(contactability))
                #print(contactability)

                #print(tensor)
                new_tensor = np.array(tensor)
                #print(new_tensor)
                for r, c in zip(new_tensor, contactability):
                    #print(r, c, end="\r")
                    r[contactability_pos] = c
                new_tensor = Tensor(new_tensor)
                new_embedding = embeddings.RelativeALEPHEmbedding.from_tensor(new_tensor, residue_embedding_class=embeddings.RelativeExpandedALEPHEmbedding0)


            else:
                new_embedding = embeddings.RelativeALEPHEmbedding.from_tensor(tensor, residue_embedding_class=embeddings.RelativeExpandedALEPHEmbedding0)
            new_embedding.name = embedding.name
            new_embedding.entity_path = embedding.entity_path
            new_embedding.sequence = embedding.sequence
            new_embedding.entity = embedding.entity
            new_embedding.param_names = embedding.param_names
            new_embedding.param_names[new_embedding.param_names.index("contactability")] = "rel_contactability"
            new_embedding.missing_indexes = embedding.missing_indexes
            #new_embedding.residue_embedding_class = embedding.residue_embedding_class
            #new_embedding.residue_embedding_name = embedding.residue_embedding_name


            new_embedding.save()
            #print(new_embedding)
            new_dataset.add(new_embedding, key=entity.name())
        except RelativeContactabilityError:
            continue

    return new_dataset











if __name__ == "__main__":
    from src.bioiain.machine import EmbeddingDataset

    dataset = EmbeddingDataset(name=f"tokens_aleph_ExpandedALEPHEmbedding0")
    dataset.load()

    calculate_relative_contactability(dataset)

