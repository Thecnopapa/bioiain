import os, json




from torch.utils.data import Dataset, DataLoader



class EmbeddingDataset(Dataset):
    def __init__(self,*args,  name, folder="./embeddings", **kwargs):
        self.name = name
        self.folder = folder
        os.makedirs(self.folder, exist_ok=True)
        self.embeddings = {}
        self.cache = None
        self.test_keys = []
        self.length = 0


    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={self.length}>"


    def __len__(self):
        return self.length

    def __getitem__(self, key):
        return self.get(key)




    def add(self, embedding, key:str|int|None=None, label_path=None):
        if key is None:
            key = len(self.embeddings)
        print("ADDING:", embedding)
        print(embedding.path)
        self.embeddings[key] = {
            "key": key,
            "start": self.length,
            "end": self.length+embedding.length,
            "embedding_path": embedding.path,
            "label_path": label_path,
            "length": embedding.length,
        }
        self.length += embedding.length
        return key



    def get(self, key, embedding=True, label=True, cache=True):

        import torch
        embedding_path = None
        label_path = None
        print("GET:", key)
        for e in self.embeddings.values():
            print(e)
            print(key < e["start"], key >= e["end"])
            if key < e["start"]: continue
            if key >= e["end"]: continue
            print("hi")
            if embedding:
                embedding_path = e["embedding_path"]
            if label:
                label_path = e["label_path"]
            break
        print("e_path", embedding_path)
        print("l_path", label_path)
        tensor = None
        label_data = None

        if self.cache is not None and cache:
            if self.cache["label_path"] == label_path:
                label_data = self.cache["label"]

            if self.cache["embedding_path"] == embedding_path:
                tensor = self.cache["tensor"]

        if tensor is None:
            if embedding_path is not None:
                tensor = torch.load(embedding_path)

        if label_data is None:
            if label_path is not None:
                if label_path.endswith(".json"):
                    label_data = json.load(open(label_path))
                elif label_path.endswith(".txt") or label_path.endswith(".label") or "." not in label_path:
                    with open(label_path, "r", encoding="utf-8") as f:
                        label_data = f.read().strip()
        print("tensor", tensor)
        print("label", label_data)
        if cache:
            self.cache = {"label_data":None, "label_path":None, "tensor":None, "embedding_path":None}
            if label_data is not None:
                self.cache["label_data"] = label_data
                self.cache["label_path"] = label_path
            if tensor is not None:
                self.cache["tensor"] = tensor.copy_
                self.cache["embedding_path"] = embedding_path

        if label and embedding: return embedding, label
        elif label: return label
        elif embedding: return tensor
        else: return None




    def add_label(self, key, label):
        self.embeddings[key]["label"] = label
        return self[key]

    def add_label_from_string(self, label, key=None):
        if key is None:
            key = len(self.embeddings) - 1
        folder = os.path.dirname(self.embeddings[key]["embedding_path"])
        fname = f"{self.name}.label.txt"

        path = os.path.join(folder, fname)

        with open(path, "w") as f:
            f.write(label)

        self.embeddings[key]["label_path"] = path
        return key

    def export(self, folder=None):
        if folder is None:
            assert self.folder is not None
            folder = self.folder
        data = {
            "name": self.name,
            "embedding_class": self[0]["embedding"].__class__.__name__,
            "list_class": self.__class__.__name__,
            "embeddings": self.embeddings,
        }
        fname = f"{self.name}.{data['list_class']}.embeddings.json"
        path = os.path.join(folder, fname)
        json.dump(data, open(path, "w"))
        return path



