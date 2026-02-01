import os, json
from operator import truediv

import torch

import torch
import torch.nn as nn

import pandas as pd

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

        os.makedirs(self.data["folder"], exist_ok=True)
        self.get_fname()

    def add_epoch(self):
        if self.data["epoch"] is None: self.data["epoch"] = 0
        else: self.data["epoch"] += 1

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

    def test(self, dataset):
        self.load()
        dataset.test()
        print(f"Testing saved model at: {self.data['path']}")
        print("Dataset:", dataset)
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

        print(f"Model EPOCH:{self.data['epoch']} correct={correct}, total={total}, accuracy={(correct / total) * 100:2.3f}%")
        #print(json.dumps(confusion, indent=4))
        df = pd.DataFrame.from_dict(confusion, orient='index')
        cf=""" Confusion Matrix: {}
            P  R  E  D  S
        T    |  {:2s}|  {:2s}|  {:2s}|  {:2s}
        R  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        U  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        T  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        H  {:2s}|{:4d}|{:4d}|{:4d}|{:4d}
        S
        """.format(
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






class MLP_MK1(CustomModel):
    def __init__(self, *args, input_dim, hidden_dims=[256 ,128], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)
        self.model = nn.Sequential(
            nn.Linear(input_dim, hidden_dims[0]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dims[0], hidden_dims[1]),
            nn.ReLU(),
            nn.Dropout(dropout),
            nn.Linear(hidden_dims[1], num_classes)
        )



    def forward(self, x):
        return self.model(x)













def model_mapping():
    mapping = {
        "internship_MLP": MLP_MK1
    }
    return mapping

