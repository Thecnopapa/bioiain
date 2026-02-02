import os, json
from copy import deepcopy


from ..utilities.logging import log
from typing import Any

from torch import Tensor
from torch.utils.data import Dataset, DataLoader


class Item(object):
    def __init__(self, tensor:Tensor, label:Any, label_to_index:dict|None=None, key:str=None, dataset=None):
        self.tensor = tensor
        self.label = label
        self.t = self.tensor
        self.l = self.label
        if label_to_index is not None:
            self.label_index = label_to_index[self.label]
            self.label_tensor = [0] * len(label_to_index)
            self.label_tensor[label_to_index[self.label]] = 1
            self.label_tensor = Tensor(self.label_tensor)
            self.li = self.label_index
            self.lt = self.label_tensor

        self.key = key
        self.dataset = dataset

    def __getitem__(self, item):
        if item in [0, "tensor", "t"]:
            return self.tensor
        elif item in [1, "label", "l"]:
            return self.label
        elif item in [2, "label_index", "li"]:
            return self.label_index
        elif item in [3, "label_tensor", "lt"]:
            return self.label_index
        else:
            raise KeyError(item)

    def __repr__(self):
        return f"Item({self.key}) T:{self.tensor.shape}, L=\"{self.label}\", from: {self.dataset}"

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i > 3:
            raise StopIteration
        return self[self.i]


class EmbeddingDataset(Dataset):
    def __init__(self,*args,  name, folder="./datasets", **kwargs):
        fname = f"{name}.dataset.json"
        path = os.path.join(folder, fname)
        self.data = dict(
            name = name,
            folder = folder,
            length = 0,
            test_length = 0,
            fname = fname,
            path = path,
            mapped = False,
        )
        self.mode="normal"
        os.makedirs(self.data["folder"], exist_ok=True)
        self.cache = None
        self.embeddings = {}
        self.splitted = {
            "test": None,
            "train": None,
        }

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.data["name"]} N={len(self)} mode={self.mode}>"


    def __len__(self):
        if self.mode == "normal": return self.data["length"]
        elif self.mode == "test": return self.splitted["test_length"]
        elif self.mode == "train": return self.splitted["train_length"]
        else: raise Exception(f"Unknown mode: {self.mode}")

    def __getitem__(self, key):
        return self.get(key)

    def __contains__(self, item):
        return item in self.embeddings.keys()

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i >= len(self):
            raise StopIteration
        r = self.get(self.i)
        self.i += 1
        return r

    def test(self):
        assert self.splitted["test"] is not None
        self.mode="test"

    def train(self):
        assert self.splitted["train"] is not None
        self.mode="train"

    def normal(self):
        self.mode="normal"

    def split(self, mode="embeddings", test_ratio=0.1, random_state=42):
        import random, math
        log(1, f"Splitting dataset...")
        if mode == "embeddings" or True:
            data = list(self.embeddings.items())

            n_keys = math.floor(len(data)*test_ratio)
            random.shuffle(data)
            test = data[0:n_keys]
            train = data[n_keys:]

            for name, dataset in zip(("test", "train"), (test, train)):
                n = 0
                self.splitted[name] = {}
                for k, s in dataset:
                    v = self.splitted[name][k] = deepcopy(s)
                    v["start"] = n
                    n += v["length"]
                    v["end"] = n
                self.splitted[name+"_length"] = n

            log(2, "(test) {} / {} (train)".format(len(self.splitted["test"]), len(self.splitted["train"])))
            return self
        else:
            raise Exception("Not implemented split method:", mode)




        # elif mode == "indices":
        #     indices = list(range(0, len(self)))
        #     n_indices = math.floor(len(self.embeddings)*test_ratio)
        #     random.shuffle(indices)
        #     indoces = indices[0:n_indices]
        #     self.test_info["indices"] = indices
        #     self.test_info["length"] = len(self.test_info["indices"])
        #     self.data["length"] -= self.test_info["length"]
        #     return self.test_info["indices"]




    def map(self) -> dict:
        log(1, "Mapping dataset...")
        self.data["label_to_index"] = {}
        self.data["index_to_label"] = {}
        for item in self:
            label = item.label
            if label in self.data["label_to_index"].keys():
                continue
            else:
                i = len(self.data["label_to_index"])
                self.data["label_to_index"][label] = i
                self.data["index_to_label"][i] = label
        self.data["mapped"] = True
        return self.data["label_to_index"]






    def add(self, embedding, key:str|int|None=None, label_path=None):
        if key is None:
            key = len(self.embeddings)
        print("ADDING:", embedding)
        #print(embedding.path)
        self.embeddings[key] = {
            "key": key,
            "start": len(self),
            "end": len(self)+embedding.length,
            "embedding_path": embedding.path,
            "label_path": label_path,
            "length": embedding.length,
            "iter_dim": embedding.iter_dim,
        }
        self.data["length"] += embedding.length
        return key



    def get(self, key, embedding=True, label=True, cache=True) -> Item:

        from torch import load as torch_load
        embedding_path = None
        label_path = None
        #print("GET:", key)
        iter_dim = 0
        rel_key = None

        if self.mode == "normal": data = self.embeddings
        elif self.mode == "test": data = self.splitted["test"]
        elif self.mode == "train": data = self.splitted["train"]
        else: raise Exception(f"Unknown mode: {self.mode}")

        if self.mode == "normal": emb_list = self.embeddings
        elif self.mode == "test": emb_list = self.splitted["test"]
        elif self.mode == "train": emb_list = self.splitted["train"]
        else: raise Exception("Not implemented split method:", self.mode)



        for e in emb_list.values():
            #print(e)
            #print(key < e["start"], key >= e["end"])
            if key < e["start"]: continue
            if key >= e["end"]: continue
            iter_dim = e["iter_dim"]

            if embedding:
                embedding_path = e["embedding_path"]
            if label:
                label_path = e["label_path"]
            rel_key = key - e["start"]
            break

        assert rel_key is not None

        #print("REL_KEY:", rel_key)

        #print("e_path", embedding_path)
        #print("l_path", label_path)
        tensor = None
        label_data = None

        if self.cache is not None and cache:
            if self.cache["label_path"] is not None:
                if self.cache["label_path"] == label_path:
                    label_data = self.cache["label_data"]

            if self.cache["embedding_path"] is not None:
                if self.cache["embedding_path"] == embedding_path:
                    tensor = self.cache["tensor"]

        if tensor is None:
            if embedding_path is not None:
                tensor = torch_load(embedding_path)

        if label_data is None:
            if label_path is not None:
                if label_path.endswith(".json"):
                    label_data = json.load(open(label_path))
                elif label_path.endswith(".txt") or label_path.endswith(".label") or "." not in label_path:
                    with open(label_path, "r", encoding="utf-8") as f:
                        label_data = f.read().strip()


        target_tensor=None
        target_label=None
        if embedding:
            target_tensor = tensor
            for i in range(iter_dim):
                target_tensor = target_tensor[0]
            target_tensor = target_tensor[rel_key]
            #print("tensor", target_tensor.shape)

        if label:
            if label_data is not None:
                target_label = label_data[rel_key]
            #print("label", target_label)

        if cache:
            self.cache = {"label_data":None, "label_path":None, "tensor":None, "embedding_path":None}
            if label_data is not None:
                self.cache["label_data"] = label_data
                self.cache["label_path"] = label_path
            if tensor is not None:
                self.cache["tensor"] = tensor
                self.cache["embedding_path"] = embedding_path

        l_to_i = None
        if "mapped" not in self.data:
            self.data["mapped"] = False # DEBUG
        if self.data["mapped"]:
            l_to_i = self.data["label_to_index"]
        return Item(target_tensor, target_label, label_to_index=l_to_i, key=key, dataset=self)




    def add_label(self, key, label):
        self.embeddings[key]["label"] = label
        return self[key]

    def add_label_from_string(self, label, key=None):
        if key is None:
            key = len(self.embeddings) - 1
        folder = os.path.dirname(self.embeddings[key]["embedding_path"])
        fname = f"{self.data['name']}.label.txt"

        path = os.path.join(folder, fname)


        with open(path, "w") as f:
            f.write(label)
        assert self.embeddings[key]["length"] == len(label)
        self.embeddings[key]["label_path"] = path
        return key

    def save(self, *args, **kwargs):
        return self.export(*args, **kwargs)

    def export(self, folder=None, save_split=False):
        if folder is None:
            assert self.data["folder"] is not None
            folder = self.data["folder"]
        data = {
            "data": self.data,
            "embeddings": self.embeddings,
        }
        if save_split:
            data["splitted"] = self.splitted
        os.makedirs(folder, exist_ok=True)
        path = os.path.join(folder, self.data["fname"])
        json.dump(data, open(path, "w"), indent=4)
        return path

    def load(self, folder=None, missing_ok=True, load_split=False):
        if folder is None:
            assert self.data["folder"] is not None
            folder = self.data["folder"]
        path = os.path.join(folder, self.data["fname"])
        if not os.path.exists(path) and missing_ok:
            log("warning", f"Dataset data not found at: {path}")
            return self
        raw = json.load(open(path, "r"))
        self.data = raw["data"]
        self.embeddings = raw["embeddings"]
        if load_split:
            try: self.splitted = raw["splitted"]
            except KeyError: log("warning", f"Dataset split info not found at: {path}")
        return self

    @classmethod
    def from_file(cls, path, load_split=False):
        raw = json.load(open(path, "r"))
        name = raw["data"]["name"]
        folder = raw["data"]["folder"]
        new = cls(name=name, folder=folder)
        new.data = raw["data"]
        new.embeddings = raw["embeddings"]
        if load_split:
            try:  new.splitted = raw["splitted"]
            except KeyError: log("warning", f"Dataset split info not found at: {path}")
        return new


class EmbeddingDataloader(DataLoader):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, collate_fn=self.collate_fn, **kwargs)


    def collate_fn(self, batch):
        pass

