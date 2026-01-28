import os, json







#BASE CLASSES

class Embedding(object):
    pass



class EmbeddingList(object):
    def __init__(self,*args,  entity=None, **kwargs):
        self.entity = entity
        self.embeddings = {}

    

    def add(self, embedding:Embedding, key:str|int|None=None, label=None):
        if key is None:
            key = str(len(self.embeddings))

        self.embeddings[key] = {
                "key": key,
                "embedding": embedding,
                "label": label,
                }

    def add_label(self, key, label):
        self.embeddings[key]["label"] = label



class ResidueEmbedding(Embedding):
    pass



class PerResidueEmbeddings(EmbeddingList):
    
    def _get_sequence(self):
        pass



class MonomerEmbedding(Embedding):
    pass



#CUSTOM CLASSES







class SaProtEmbeddings(PerResidueEmbeddings):
    
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.sequence = None
        self.fs_tokens = None


    def generate_embeddings(self, *args, **kwargs):


        self._get_sequence()
        self._get_fs()

        self._run_saprot(self.sequence, self.fs_tokens)


    def _run_saprot(self, sequence, fs_tokens):
        pass

    def _get_fs(self):
        pass
        






class MonomerInterfaceEmbeddings(EmbeddingList):

    def __init__(self, monomer):
        super().__init__(self, monomer)

    
    def saprot_ebeddings(self):
        for residue in self.entity.get_residues():
            print(residue)





    






