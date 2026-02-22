import os, json, platform

from sklearn.metrics import confusion_matrix

from ..utilities.logging import *
from ..utilities.parallel import *

import torch
import torch.nn as nn

import pandas as pd
import numpy as np

from .datasets import Item, Dataset

from torch.utils.tensorboard import SummaryWriter
import datetime
import psutil


class ModelNotFound(Exception):
    pass




class CustomLoss(object):
    pass

class CustomModel(nn.Module):
    def __init__(self, name, in_shape, lr=0.001, folder="./models"):
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
        self.data["in_shape"] = in_shape
        os.makedirs(self.data["folder"], exist_ok=True)

        self.criterions = {
            "default": nn.MSELoss()
        }
        self.optimisers = {
            "default": {
                "class":torch.optim.Adam,
                "kwargs":{"lr":lr},
            }
        }

        self.layers = {
            "default": {
            }
        }
        self.submodels = {
        }
        self.running_loss = {"total":0, "default":0}

        log("header", f"Model initialised: {self}")




    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.data['name']} loss={self.criterions['default'].__class__.__name__} mode={self.mode} epoch={self.data['epoch']}>"


    def _mount_submodels(self):
        log(1, "Mounting submodels...")
        for k, layer_set in self.layers.items():
            self.submodels[k] = nn.Sequential(*layer_set.values())
            self.optimisers[k] = self.optimisers[k]["class"](self.submodels[k].parameters(), **self.optimisers[k]["kwargs"])
        self.mounted = True

        self.reset_loss()

        self.writer = SummaryWriter(log_dir=f"runs/{self.data['name']}/{self.__class__.__name__}/{self.optimisers["default"].__class__.__name__}-{self.criterions["default"].__class__.__name__}-{datetime.datetime.now()}")
        #self.writer.add_graph(self, torch.rand(self.data["in_shape"]))


        return self.submodels.keys()

    def reset_loss(self):
        self.running_loss["total"] = 0
        for c in self.running_loss.keys():
            self.running_loss[c] = 0

    def send_run(self, *args, **kwargs):
        if self.writer is not None:
            folder = platform.node()
            run = writer.name
            file = writer.path
            return send_tensorboard_run(*args, folder=folder, run=run, file=file, **kwargs)



    def write_loss(self):
        for c, rl in self.running_loss.items():
            total = self.running_loss["total"]
            if isinstance(rl, torch.Tensor):
                rl = rl.item()
            #print(total, c, rl)

            if total == 0: av_loss = rl
            else: av_loss = rl / total
            if self.writer is not None:
                #print("writing", av_loss)
                self.writer.add_scalar(f"loss/{c}", float(av_loss), self.data["epoch"])


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

        self.write_loss()
        self.reset_loss()
        if self.data["epoch"] is None: self.data["epoch"] = 1
        else: self.data["epoch"] += 1
        return self.data["epoch"]

    def get_fname(self, add_epoch=False) -> str:
        if add_epoch and self.data["epoch"] is not None:
            fname = f"{self.data['model']}_{self.data['name']}_E{self.data['epoch']}"
        else:
            fname = f"{self.data['model']}_{self.data['name']}"
        return fname

    def save(self, path=None, add_epoch=False, temp=False):
        if path is None:
            path = os.path.join(self.data["folder"], self.get_fname(add_epoch=add_epoch))
        if not path.endswith(".model.pt"):
            model_path = path + ".model.pt"
        if temp:
            model_path = model_path.replace(".model.pt", ".temp.model.pt")
        self.data["path"] = model_path
        torch.save(self.state_dict(), model_path)
        return self.export(path=model_path, add_epoch=add_epoch)

    def export(self, path=None, add_epoch=False):
        if path is None:
            path = os.path.join(self.data["folder"], self.get_fname(add_epoch=add_epoch))
        if path.endswith(".model.pt"):
            data_path = path.replace(".model.pt", ".data.json")
        elif not path.endswith(".data.json"):
            data_path = path + ".data.json"
        else:
            data_path = path
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
        if dataset.data.get("lab_count", False):
            self.data["lab_count"] = dataset.data["lab_count"]

        return self.data["label_to_index"], dataset.data["index_to_label"]

    def test(self, dataset, re_load=True):
        if re_load:
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
            truths = []
            preds = []

            confusion = {k: {l:0 for l in label_to_index.keys()} for k in label_to_index.keys()}

            for n, item in enumerate(dataset):
                if not is_cluster:
                    print(f"\033]0;Testing {(total/len(dataset))*100:3.0f}%\a", end="\r")
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
                        preds.append(p)
                        truths.append(l)

                        if not is_cluster:
                            print(f"{n:4d}/{len(dataset):4d} PRED: {pred.item()}, TRUTH: {truth}, CORRECT: {pred.item() == truth}", end="\r")
                        if pred == truth:
                            correct += 1
                else:
                    if abs(l-out.item()) <= 0.05:
                        correct += 1

        print()
        print(json.dumps(confusion, indent=4))
        try:
            self.writer.add_scalar(f"accuracy/dual/total", (correct / total) * 100, self.data["epoch"])
            self.writer.add_scalar(f"accuracy/dual/contactability", (confusion["contactability"]["right"] / total) * 100, self.data["epoch"])
            self.writer.add_scalar(f"accuracy/dual/outer", (confusion["outer"]["right"] / total) * 100, self.data["epoch"])
        except:
            pass


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
        if len(preds) > 0 and len(truths) > 0:
            from ..visualisation.plots import plot_confusion
            plot_confusion(preds, truths, title=f"{self}", classes = label_to_index.keys())
        self.save(temp=not re_load)


    def loss(self,
             output:torch.Tensor,
             item:Item,
             criterion_name:str="mode",
             backwards:bool=True,
             zero_optims:str|None="mode") -> torch.Tensor|float:
        self.zero_grad(zero_optims)

        if criterion_name == "mode": criterions = [self.mode]

        elif criterion_name == "all":
            criterions = [n for n in self.criterions.keys()]
        else:
            criterions = [criterion_name]


        losses = []

        for criterion in criterions:
            if isinstance(self.criterions[criterion], CustomLoss):
                losses.append(self.criterions[criterion](output, item))

            elif hasattr(item, "lt"):
                #print("LT", item.lt)
                losses.append(self.criterions[criterion](output, item.lt))
            else:
                #print("L", item.l)
                try:
                    losses.append(self.criterions[criterion](output, torch.Tensor([item.l])))
                except:
                    print(item)
                    print(item.l)
                    print(item.__dict__)
                    raise

        if len(losses) > 1:
            loss = torch.sum(losses)
        else:
            loss = losses[0]


        if backwards:
            self.running_loss["total"] += 1
            for l,c in zip(losses, criterions):
                self.running_loss[c] += l

            loss.backward()


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








class DUAL_MLP_MK5(CustomModel):
    def __init__(self, *args, hidden_dims=[2560, 1280, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], hidden_dims[0]),
            "relu3": nn.LeakyReLU(),
            "l4": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "relu4": nn.LeakyReLU(),
            "l5": nn.Linear(hidden_dims[1], hidden_dims[0]),
            "relu5": nn.LeakyReLU(),
            "l6": nn.Linear(hidden_dims[0], hidden_dims[-1]),
            "drop2": nn.Dropout(dropout),
            "last": nn.Linear(hidden_dims[-1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }



        self.criterions["default"] = self.CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])

        self._mount_submodels()



    class CustomHalfHalf(CustomLoss):
        def __init__(self, weights=None):
            log(1, "Using label weights:")

            if type(weights) is dict:
                weights = weights.values()
            w = np.array(list(weights))
            w = w
            w = w / w.sum()
            w = torch.Tensor(w)
            log(2, w)
            self.CEL = nn.MSELoss()
            self.weight = w

        def __call__(self, o, item):
            true_index = item.li
            true_tensor = item.lt
            pred = torch.max(o, dim=0)[1]
            if item.li < 5:
                true_tensor[:5] = 0.5 
            else:
                true_tensor[5:] = 0.5 
            true_tensor[item.li:item.li+1] = 1.

            #print(true_tensor)
            weighted_out = o# * self.weight
            #print("Weighted out:", weighted_out)
            loss = self.CEL(o, true_tensor)
            #if true_index < 5 == torch.max(o, dim=0)[1] < 5:
            #    loss *= 0.5
            loss *= self.weight[pred]
            return loss










class DUAL_MLP_MK4(CustomModel):
    def __init__(self, *args, hidden_dims=[2560, 128], num_classes=4, dropout=0.2, weights=None, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.LeakyReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.LeakyReLU(),
            "l3": nn.Linear(hidden_dims[1], num_classes),
            "softmax": nn.Softmax(dim=0)
        }






        self.criterions["default"] = self.CustomHalfHalf(weights)
        self.data["weights"] = list([w.item() for w in self.criterions["default"].weight])

        self._mount_submodels()



    class CustomHalfHalf(CustomLoss):
        def __init__(self, weights=None):
            log(1, "Using label weights:")

            if type(weights) is dict:
                weights = weights.values()
            w = np.array(list(weights))
            w = w
            w = w / w.sum()
            w = torch.Tensor(w)
            log(2, w)
            self.CEL = nn.MSELoss()
            self.weight = w

        def __call__(self, o, item):
            true_index = item.li
            true_tensor = item.lt
            pred = torch.max(o, dim=0)[1]
            if item.li < 5:
                true_tensor[:5] = 0.5 
            else:
                true_tensor[5:] = 0.5 
            true_tensor[item.li:item.li+1] = 1.

            #print(true_tensor)
            weighted_out = o# * self.weight
            #print("Weighted out:", weighted_out)
            loss = self.CEL(o, true_tensor)
            #if true_index < 5 == torch.max(o, dim=0)[1] < 5:
            #    loss *= 0.5
            loss *= self.weight[pred]
            return loss





class DUAL_MLP_MK3(CustomModel):
    def __init__(self, *args, hidden_dims=[640, 1280, 128], num_classes=4, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)

        self.data["num_classes"] = num_classes
        self.data["hidden_dims"] = hidden_dims
        self.data["dropout"] = dropout

        self.layers["default"] = {
            "l1": nn.Linear(self.data["in_shape"][0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.ReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.ReLU(),
            "l3": nn.Linear(hidden_dims[1], hidden_dims[2]),
            "drop3": nn.Dropout(dropout),
            "relu3": nn.ReLU(),
            "l4": nn.Linear(hidden_dims[2], num_classes),
            "softmax": nn.Softmax(dim=0)
        }

        #self.criterions["default"] = self.DualLoss(self)
        self.criterions["default"] = nn.CrossEntropyLoss()

        self._mount_submodels()






class DUAL_MLP_MK2(CustomModel):
    def __init__(self, *args, hidden_dims=[2560, 1280, 128], num_classes=2, dropout=0.2, **kwargs):
        super().__init__(*args, **kwargs)



        self.layers["default"] = {
            "l1": nn.Linear(self.in_shape[0], hidden_dims[0]),
            "drop1": nn.Dropout(dropout),
            "relu1": nn.ReLU(),
            "l2": nn.Linear(hidden_dims[0], hidden_dims[1]),
            "drop2": nn.Dropout(dropout),
            "relu2": nn.ReLU(),
            "l3": nn.Linear(hidden_dims[1], hidden_dims[2]),
            "drop3": nn.Dropout(dropout),
            "relu3": nn.ReLU(),
            "l4": nn.Linear(hidden_dims[2], num_classes),
            "hardtahn": nn.Hardtanh(min_val=0.)
        }

        self.criterions["default"] = self.DualLoss(self)

        self._mount_submodels()


    class DualLoss(object):
        def __init__(self, model):
            self.writer = model.writer
            self.model = model

        def __name__(self):
            return "DualLoss"

        def __call__(self, o, t):
            true_contact, true_outer = t[0], t[1]
            out_contact, out_outer = o[0], o[1]

            outer_loss = abs(true_outer - out_outer)
            contact_loss = abs(true_contact - out_contact)
            if "outer" not in self.model.running_loss: self.model.running_loss["outer"] = 0
            if "contactability" not in self.model.running_loss: self.model.running_loss["contactability"] = 0
            self.model.running_loss["outer"] += outer_loss
            self.model.running_loss["contactability"] += contact_loss
 
            return contact_loss * outer_loss


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

        self.criterions["default"] = self.DualLoss(self)

        self._mount_submodels()


    class DualLoss(object):
        def __init__(self, model):
            self.writer = model.writer
            self.model = model

        def __name__(self):
            return "DualLoss"

        def __call__(self, o, t):
            true_contact, true_outer = t[0], t[1]
            out_contact, out_outer = o[0], o[1]

            outer_loss = abs(true_outer - out_outer)
            contact_loss = abs(true_contact - out_contact)
            if "outer" not in self.model.running_loss: self.model.running_loss["outer"] = 0
            if "contactability" not in self.model.running_loss: self.model.running_loss["contactability"] = 0
            self.model.running_loss["outer"] += outer_loss
            self.model.running_loss["contactability"] += contact_loss
 
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


