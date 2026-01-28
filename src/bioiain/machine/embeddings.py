import os, json







#BASE CLASSES

class Embedding(object):

    def __init__(self, *args, name=None, folder=None, **kwargs):
        if folder is None:
            fodler = "./embeddings"
        self.folder =folder
        self.name=name
        self.path = None

    def from_file(path):
        self.name = path.split(".")[0]
        self.path = path
        self.folder = os.path.dirname(path)





class EmbeddingList(object):
    def __init__(self,*args,  name, folder=None, **kwargs):
        self.name = name
        self.folder = folder
        
        self.embeddings = {}

    

    def add(self, embedding:Embedding, key:str|int|None=None, label=None):
        if key is None:
            key = str(len(self.embeddings))

        self.embeddings[key] = {
                "key": key,
                "embedding_path": embedding.path,
                "label": label,
                }
        return self[key]

    def __getitem__(self, key):
        return self.embeddings[key]

    def add_label(self, key, label):
        self.embeddings[key]["label"] = label
        return self[key]

    def export(self, folder=None):
        if folder is None:
            assert self.folder is not None
            folder = self.folder
        data = {
                "name": self.name,
                "embeddint_class": self[0]["embedding"].__class__.__name__,
                "list_class": self.__class__.__name__,
                "embeddings" = self.embeddings,
                }
        fname = f"{self.name}.{data['list_class']}.embeddings.json"
        path = os.path.join(fodler, fname)
        json.dump(data, open(path, "w"))
        return path



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





    






