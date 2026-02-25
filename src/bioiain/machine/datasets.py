import os, json
from copy import deepcopy


from ..utilities.logging import log
from ..utilities.exceptions import *
from typing import Any

from torch import Tensor
from torch.utils.data import Dataset, DataLoader


class Item(object):
    def __init__(self, tensor:Tensor, label:Any, label_to_index:dict|None=None, key:str=None, dataset=None):
        self.tensor = tensor
        self.label = label
        self.t = self.tensor
        self.l = self.label
        if label_to_index is not None and len(label_to_index) > 1:
            if type(self.label) in [int, str]:
                #print("LABEL IS INT/STR")
                self.label_index = label_to_index[self.label]
                self.label_tensor = [0] * len(label_to_index)
                self.label_tensor[label_to_index[self.label]] = 1
                self.label_tensor = Tensor(self.label_tensor)
                self.li = self.label_index
                self.lt = self.label_tensor
            elif type(self.label) in (list, tuple):
                #print("LABEL IS LIST/TUPLE")
                self.label_tensor = Tensor(self.label)
                self.lt = self.label_tensor
        #print(f"LABEL IS {type(self.label)}", type(self.label) in (list, tuple), label_to_index is not None , len(label_to_index) > 1)


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
        return f"<bi.{self.__class__.__name__}:{self.key} T:{self.tensor.shape}, L=\"{self.label}\", from: {self.dataset}>"

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i > 3:
            self.i = None
            raise StopIteration
        return self[self.i]

    def to(self, device):
        try:
            self.l = self.l.to(device)
        except:
            pass
        self.t = self.t.to(device)
        if hasattr(self, "lt"):
            self.lt = self.lt.to(device)
        if hasattr(self, "label_tensor"):
            self.label_tensor = self.label_tensor.to(device)


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
            label_key = "label_path",
            deleted_indexes = 0,
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
        if self.data["deleted_indexes"] > 0:
            return f"<bi.{self.__class__.__name__}:{self.data["name"]} N={len(self)} mode={self.mode} deleted={self.data.get('deleted',False)}>"
        else:
            return f"<bi.{self.__class__.__name__}:{self.data["name"]} N={len(self)} mode={self.mode}>"


    def __len__(self):
        if self.mode == "normal": 
            return self.data["length"] - self.data["deleted_indexes"]
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
        try:
            r = self.get(self.i)
        except DeletedIndex:
            self.i += 1
            r = self.__next__()
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


    def use_label(self, label_key):
        self.data["label_key"] = label_key


    def split(self, mode="embeddings", test_ratio=0.1, random_state=42):
        import random, math
        log(1, f"Splitting dataset...")
        if mode == "embeddings" or True:
            data = [e for e in self.embeddings.items() if not e[1].get("deleted", False)]

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


    def map(self, single_lab=False, label_to_index:dict|None=None, reuse=True) -> dict:
        log(1, "Mapping dataset...")
        
        if self.data["mapped"] and (self.data.get("mapped_label", None) == self.data["label_key"]) and reuse:
            return self.data["label_to_index"]

        self.data["mapped"] = False
        self.data["mapped_label"] = None
        self.data["label_to_index"] = {}
        self.data["index_to_label"] = {}
        self.data["lab_count"] = {}

        lab_count = {}


        if label_to_index is not None:
            labels = label_to_index.keys()

        elif single_lab:
            labels = [0]

        else:
            labels = []
            for item in self:
                label = item.label
                #print(label, item)
                if label not in labels:
                    labels.append(label)
                    lab_count[label] = 0
                lab_count[label] += 1

        for n, k in enumerate(sorted(labels)):
                self.data["label_to_index"][str(k)] = int(n)
                self.data["index_to_label"][int(n)] = str(k)
                if len(lab_count) > 0:
                    self.data["lab_count"][str(k)] = int(lab_count[k])

        self.data["mapped_label"] = self.data["label_key"]
        self.data["mapped"] = True
        return self.data["label_to_index"]


    def add(self, embedding, key:str|int|None=None, label_path=None, fasta=False):
        if key is None:
            key = len(self.embeddings)
        #print("ADDING:", embedding)
        #print(embedding.path)
        self.embeddings[key] = {
            "key": key,
            "start": len(self),
            "end": len(self)+embedding.length,
            "embedding_path": embedding.path,
            "label_path": label_path,
            "length": embedding.length,
            "iter_dim": embedding.iter_dim,
            "deleted": False,

        }
        if embedding.keep_padding:
            self.embeddings[key]["padding"]= embedding.padding,
        else:
            self.embeddings[key]["padding"] = 0
        self.data["length"] += embedding.length
        if fasta and hasattr(embedding, "sequence"):
            self._add_to_fasta(key, embedding.sequence)

        self.data["mapped"] = False
        return key


    def remove(self, key):
        self.data["deleted_indexes"] += self.embeddings[key]["end"] - self.embeddings[key]["start"]
        self.embeddings[key]["deleted"] = True
        self.data["mapped"] = False


    def _add_to_fasta(self, key, sequence):
        if "fasta_path" in self.data:
            fasta_path = self.data["fasta_path"]
            mode = "a"
        else:
            fasta_path = os.path.join(self.data["folder"], self.data["name"]+".dataset.fasta")
            self.data["fasta_path"] = fasta_path
            mode = "w"
        with open(fasta_path, mode) as f:
            if mode == "w":
                #f.write(f"# FASTA for dataset: {self.data["name"]}\n\n")
                pass
            f.write(f"> {key}\n")
            f.write(f"{sequence}\n\n")


    def get(self, key, embedding=True, label=True, cache=True, label_key=None) -> Item:
        from torch import load as torch_load

        if label_key is None:
            label_key = self.data["label_key"]
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

            if e.get("deleted", False):
                raise DeletedIndex(f"Index: {key}")

            iter_dim = e["iter_dim"]

            if embedding:
                embedding_path = e["embedding_path"]
            if label:
                label_path = e[label_key]
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
                elif label_path.endswith(".csv"):
                    with open(label_path, "r") as f:
                        label_data = [l for l in f.read().strip().split(",")]
                        for n, l in enumerate(label_data):
                            if ":" in l:
                                label_data[n] = [float(ll) for ll in l.split(":")]
                            else:
                                try:
                                    label_data[n] = float(l)
                                except ValueError:
                                    label_data[n] = l.strip()

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
        #print(l_to_i)
        return Item(target_tensor, target_label, label_to_index=l_to_i, key=key, dataset=self)


    def add_label(self, key, label):
        self.embeddings[key]["label"] = label
        self.data["mapped"] = False

        return self[key]


    def add_label_from_string(self, label, key=None, var_name="label_path"):
        if key is None:
            key = len(self.embeddings) - 1

        if self.embeddings[key]["padding"] > 0:
            label = ">"*self.embeddings[key]["padding"] + label + "<"*self.embeddings[key]["padding"]
        folder = os.path.dirname(self.embeddings[key]["embedding_path"])
        fname = f"{var_name}.label.txt"

        path = os.path.join(folder, fname)


        with open(path, "w") as f:
            f.write(label)
        assert self.embeddings[key]["length"] == len(label)
        self.embeddings[key][var_name] = path
        self.data["mapped"] = False

        return key


    def add_label_from_list(self, label, key=None, var_name="label_path", paddings=(None, None)):
        if key is None:
            key = len(self.embeddings) - 1
        folder = os.path.dirname(self.embeddings[key]["embedding_path"])
        fname = f"{var_name}.label.csv"

        path = os.path.join(folder, fname)

        labels = []
        for l in label:
            if type(l) in [list, tuple]:
                labels.append(":".join([str(ll) for ll in l]))
            else:
                labels.append(str(l))


        if self.embeddings[key]["padding"] > 0:
            label = paddings[0]*self.embeddings[key]["padding"] + label + paddings[1]*self.embeddings[key]["padding"]


        with open(path, "w") as f:
            f.write(",".join(labels))
        assert self.embeddings[key]["length"] == len(label)
        self.embeddings[key][var_name] = path
        self.data["mapped"] = False
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

