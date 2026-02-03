import os, json
from ..utilities.logging import log

import torch
import torch.nn as nn

import pandas as pd

from .datasets import Item, Dataset


class ModelNotFound(Exception):
    pass

class CustomModel(nn.Module):
    def __init__(self, name, folder="./models"):
        super().__init__()
        self.data = {}
        self.data["name"] = name
        self.data["folder"] = folder
        self.data["epoch"] = None
        self.data["path"] = False
        self.data["model"] = self.__class__.__name__
        self.mode = "default"
        self.mounted = False
        os.makedirs(self.data["folder"], exist_ok=True)

        self.criterions = {
            "default": nn.MSELoss()
        }
        self.optimizers = {
            "default": {
                "class":torch.optim.Adam,
                "kwargs":{"lr":0.001},
            }
        }

        self.layers = {
            "default": {
            }
        }
        self.submodels = {
        }
        log("header", f"Model initialised: {self}")




    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.data['name']} mode={self.mode} epoch={self.data['epoch']}>"


    def _mount_submodels(self):
        log(1, "Mounting submodels...")
        for k, layer_set in self.layers.items():
            self.submodels[k] = nn.Sequential(*layer_set.values())
            self.optimizers[k] = self.optimizers[k]["class"](self.submodels[k].parameters(), **self.optimizers[k]["kwargs"])
        self.mounted = True
        return self.submodels.keys()


    def set_epoch(self, epoch):
        self.data["epoch"] = epoch
        return self.data["epoch"]

    def add_epoch(self):
        if self.data["epoch"] is None: self.data["epoch"] = 1
        else: self.data["epoch"] += 1
        return self.data["epoch"]

    def get_fname(self, add_epoch=False) -> str:
        if add_epoch and self.data["epoch"] is not None:
            fname = f"{self.data['model']}_{self.data['name']}_E{self.data['epoch']}"
        else:
            fname = f"{self.data['model']}_{self.data['name']}"
        return fname

    def save(self, path=None, add_epoch=False):
        if path is None:
            path = os.path.join(self.data["folder"], self.get_fname(add_epoch=add_epoch))
        model_path = path + ".model.pt"
        self.data["path"] = model_path
        torch.save(self.state_dict(), model_path)
        return self.export(path=path, add_epoch=add_epoch)

    def export(self, path=None, add_epoch=False):
        if path is None:
            path = os.path.join(self.data["folder"], self.get_fname(add_epoch=add_epoch))
        data_path = path + ".data.json"
        json.dump(self.data, open(data_path, "w"), indent=4)
        return data_path

    def load(self, data_path=None, epoch=None, weights_only=False):
        if data_path is None:
            data_path = os.path.join(self.data["folder"], self.get_fname(add_epoch=epoch))+".data.json"
        if not os.path.exists(data_path):
            raise ModelNotFound(data_path)
        raw_data = json.load(open(data_path, "r"))
        self.data = raw_data
        self.load_state_dict(torch.load(self.data["path"], weights_only=weights_only))
        return self

    def add_map(self, dataset):
        self.data["label_to_index"] = dataset.data["label_to_index"]
        self.data["index_to_label"] = dataset.data["index_to_label"]

        return self.data["label_to_index"], dataset.data["index_to_label"]

    def test(self, dataset):
        self.load()
        dataset.test()
        log(1, f"Testing: using model saved at: {self.data['path']}")
        log(2, "Dataset:", dataset)
        #print(json.dumps(dataset.splitted["test"], indent=4))
        # for e in dataset.splitted["test"].values():
        #     print(e)
        #     print()
        #     label = ""
        #     for n in range(e["length"]):
        #         label += dataset[n][1]
        #     print("\n"+label)
        #     exit()

        label_to_index = dataset.data["label_to_index"]
        index_to_label = dataset.data["index_to_label"]
        with torch.no_grad():
            correct = 0
            total = 0

            confusion = {k: {l:0 for l in label_to_index.keys()} for k in label_to_index.keys()}

            for item in dataset:
                # truth = [0] * len(label_to_index)
                # truth[label_to_index[item.l]] = 1
                l = item.l
                truth = label_to_index[l]



                out = self(item.t)
                pred = out.argmax(dim=0)
                p = index_to_label[pred.item()]
                confusion[l][p] += 1

                print(f"PRED: {pred.item()}, TRUTH: {truth}, CORRECT: {pred.item() == truth}", end="\r")
                total += 1
                if pred == truth:
                    correct += 1

        log(1, f"Results: EPOCH:{self.data['epoch']} correct={correct}, total={total}, accuracy={(correct / total) * 100:2.3f}%")
        #print(json.dumps(confusion, indent=4))
        df = pd.DataFrame.from_dict(confusion, orient='index')
        cf=""" ** Confusion Matrix: {} **
            P  R  E  D  S
        T    |  {:2s}|  {:2s}|  {:2s}|  {:2s}
        R  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        U  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        T  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        H  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        S
********""".format(
            f"EPOCH:{self.data['epoch']} correct={correct}, total={total}, accuracy={(correct / total) * 100:2.3f}%",
            *confusion.keys(),
            list(confusion.keys())[0], *confusion[list(confusion.keys())[0]].values(),
            list(confusion.keys())[1], *confusion[list(confusion.keys())[1]].values(),
            list(confusion.keys())[2], *confusion[list(confusion.keys())[2]].values(),
            list(confusion.keys())[3], *confusion[list(confusion.keys())[3]].values(),
        )
        print(cf)
        #df.rename({0:"Truth\\Pred"}, inplace = True)
        #print(df)
        self.data["confusion_matrix"] = cf
        self.save()


    def loss(self,
             output:torch.Tensor,
             item:Item,
             criterion_name:str="mode",
             backwards:bool=True,
             zero_optims:str|None="mode") -> torch.Tensor|float:

        self.zero_grad(zero_optims)

        if criterion_name == "mode": criterion_name = self.mode

        if criterion_name == "all":
            loss = torch.sum([self.criterions[criterion_name](output, item.lt) for criterion_name in self.criterions.keys()])
        else:
            loss = self.criterions[criterion_name](output, item.lt)
        if backwards:
            loss.backward()
        return loss

    def step(self, optimizer_name:str|None="mode") -> bool:
        if optimizer_name is None: return False
        if optimizer_name == "mode": optimizer_name = self.mode

        if optimizer_name == "all":
            for optimizer in self.optimizers.values():
                optimizer.step()
        else:
            self.optimizers[optimizer_name].step()
        return True

    def zero_grad(self, optimizer_name:str|None="mode") -> bool:
        if optimizer_name is None: return False
        if optimizer_name == "mode": optimizer_name = self.mode

        if optimizer_name == "all":
            for optimizer in self.optimizers.values():
                optimizer.zero_grad()
        else:
            self.optimizers[optimizer_name].zero_grad()
        return True

    def set_mode(self, mode:str):
        self.mode = mode
        return self.mode

    def __call__(self, x, **kwargs):
        return self.forward(x, **kwargs)


    def forward(self, x, submodel_name="mode"):
        if submodel_name is None: return None
        if submodel_name == "mode": submodel_name = self.mode

        if submodel_name == "all":
            for submodel in self.submodels.values():
                x = submodel(x)
        else:
            x = self.submodels[submodel_name](x)
        return x


class MLP_MK2(CustomModel):
    def __init__(self, *args, input_dim, hidden_dims=[256 ,128], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(input_dim, hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            #"drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            #"drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }

        self._mount_submodels()


class MLP_MK1(CustomModel):
    def __init__(self, *args, input_dim, hidden_dims=[256 ,128], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(input_dim, hidden_dims[0]),
            "relu1": nn.ReLU(),
            "drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.ReLU(),
            "drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes)
        }

        self._mount_submodels()


