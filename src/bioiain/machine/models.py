import os, json

from sklearn.metrics import confusion_matrix

from ..utilities.logging import log

import torch
import torch.nn as nn

import pandas as pd

from .datasets import Item, Dataset

from torch.utils.tensorboard import SummaryWriter
import datetime


class ModelNotFound(Exception):
    pass

class CustomModel(nn.Module):
    def __init__(self, name, in_shape, folder="./models"):
        super().__init__()
        self.data = {}
        self.data["name"] = name
        self.data["folder"] = folder
        self.data["epoch"] = None
        self.data["path"] = False
        self.data["model"] = self.__class__.__name__
        self.mode = "default"
        self.writer = None
        self.mounted = False
        self.in_shape = in_shape
        os.makedirs(self.data["folder"], exist_ok=True)

        self.criterions = {
            "default": nn.MSELoss()
        }
        self.optimisers = {
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
            self.optimisers[k] = self.optimisers[k]["class"](self.submodels[k].parameters(), **self.optimisers[k]["kwargs"])
        self.mounted = True

        self.writer = SummaryWriter(log_dir=f"runs/{self.data['name']}/{self.optimisers["default"].__class__.__name__}-{self.criterions["default"].__class__.__name__}-{datetime.datetime.now()}")
        #self.writer.add_graph(self, torch.rand(self.in_shape))


        return self.submodels.keys()


    def set_epoch(self, epoch):
        self.data["epoch"] = epoch
        return self.data["epoch"]

    def add_epoch(self):

        if self.writer is not None:
            for set_name, layers in self.layers.items():
                for layer_name, layer in layers.items():
                    if hasattr(layer, "weight"):
                        self.writer.add_histogram(f"{set_name}/weight/{layer_name}", layer.weight, self.data["epoch"])
                    if hasattr(layer, "bias"):
                        self.writer.add_histogram(f"{set_name}/bias/{layer_name}", layer.bias,  self.data["epoch"])

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
        with (torch.no_grad()):
            correct = 0
            total = 0

            confusion = {k: {l:0 for l in label_to_index.keys()} for k in label_to_index.keys()}

            for n, item in enumerate(dataset):
                l = item.l

                if n == 0:
                    if type(l) in (list, tuple):
                        confusion = {
                            k: {"right": 0, "wrong": 0} for k in label_to_index.keys()
                        }

                out = self(item.t)
                total += 1
                if len(label_to_index) > 1:
                    if type(l) in (list, tuple):
                        #print(l)
                        truth_contact, truth_outer = l
                        truth_outer = truth_outer > 0.5
                        out_contact, out_outer = out
                        out_contact, out_outer = out_contact.item(), out_outer.item() > 0.5
                        #print("OUTER", truth_outer, out_outer, truth_outer == out_outer)
                        if n == 0:
                            confusion["outer"]["pred_inner"] = 0
                            confusion["outer"]["pred_outer"] = 0
                        if out_outer:
                            confusion["outer"]["pred_outer"] += 1
                        else:
                            confusion["outer"]["pred_inner"] += 1

                        if truth_outer == out_outer:
                            confusion["outer"]["right"] +=1
                            correct += 0.5
                        else:
                            confusion["outer"]["wrong"] += 1


                        #print("CONTACTABILITY", truth_contact, out_contact, abs(truth_contact - out_contact) <= 0.1)

                        if abs(truth_contact - out_contact) <= 0.1:
                            confusion["contactability"]["right"] += 1
                            correct += 0.5
                        else:
                            confusion["contactability"]["wrong"] += 1
                        #print(confusion)


                    else:
                        truth = label_to_index[l]
                        pred = out.argmax(dim=0)
                        p = index_to_label[pred.item()]
                        confusion[l][p] += 1

                        print(f"PRED: {pred.item()}, TRUTH: {truth}, CORRECT: {pred.item() == truth}", end="\r")
                        if pred == truth:
                            correct += 1
                else:
                    if abs(l-out.item()) <= 0.05:
                        correct += 1


        print(json.dumps(confusion, indent=4))
        log(1, f"Results: EPOCH:{self.data['epoch']-1} correct={correct}, total={total}, accuracy={(correct / total) * 100:2.3f}%")
        #print(json.dumps(confusion, indent=4))
        if len(label_to_index) == 4:
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
            #print(item, item.lt)
            if hasattr(item, "lt"):
                #print("LT", item.lt)
                loss = self.criterions[criterion_name](output, item.lt)
            else:
                #print("L", item.l)
                loss = self.criterions[criterion_name](output, torch.Tensor([item.l]))
        if backwards:
            loss.backward()

        if self.writer is not None:
            self.writer.add_scalar(f"loss/{criterion_name}", loss, self.data["epoch"])


        return loss

    def step(self, optimizer_name:str|None="mode") -> bool:
        if optimizer_name is None: return False
        if optimizer_name == "mode": optimizer_name = self.mode

        if optimizer_name == "all":
            for optimizer in self.optimisers.values():
                optimizer.step()
        else:
            self.optimisers[optimizer_name].step()
        return True

    def zero_grad(self, optimizer_name:str|None="mode") -> bool:
        if optimizer_name is None: return False
        if optimizer_name == "mode": optimizer_name = self.mode

        if optimizer_name == "all":
            for optimizer in self.optimisers.values():
                optimizer.zero_grad()
        else:
            self.optimisers[optimizer_name].zero_grad()
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


class DUAL_MLP_MK1(CustomModel):
    def __init__(self, *args, hidden_dims=[128, 256], num_classes=2, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            #"drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            #"drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            # "softmax": nn.Softmax(dim=0)
        }

        self.criterions["default"] = self.dual_loss

        self._mount_submodels()

    def dual_loss(self, o, t):
        #print("CALCULATING DUAL LOSS")
        #print("OUT",o)
        #print("TRUTH",t)
        true_contact, true_outer = t[0], t[1]
        out_contact, out_outer = o[0], o[1]


        outer_loss = abs(true_outer - out_outer)
        contact_loss = abs(true_contact - out_contact)
        #print("LOSSES")
        #print(contact_loss, outer_loss)
        self.writer.add_scalar(f"loss/dual/contact", contact_loss, self.data["epoch"])
        self.writer.add_scalar(f"loss/dual/outer", outer_loss, self.data["epoch"])
        return contact_loss + outer_loss


class MLP_MK3(CustomModel):
    def __init__(self, *args, hidden_dims=[128 ,256], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            "drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            "drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            #"softmax": nn.Softmax(dim=0)
        }


        self.criterions["default"] = self.simpleloss

        self._mount_submodels()

    @staticmethod
    def simpleloss(o, t):
            return abs(o-t)

class MLP_MK2(CustomModel):
    def __init__(self, *args, hidden_dims=[256 ,128], num_classes=8, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "relu1": nn.LeakyReLU(),
            #"drop1": nn.Dropout(dropout),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            #"drop2": nn.Dropout(dropout),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            #"softmax": nn.Softmax(dim=0)
        }


        self.criterions["default"] = self.simpleloss

        self._mount_submodels()

    @staticmethod
    def simpleloss(o, t):
            return abs(o-t)

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


