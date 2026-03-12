import os, json, platform

import torchvision.transforms.v2.functional
from sklearn.metrics import confusion_matrix
import PIL

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
from . import DEVICE

from ..utilities.exceptions import *
from .losses import *



class BaseModel(nn.Module):
    def __init__(
            self, 
            name:str, 
            in_shape:list|tuple, 
            lr:float=0.001, 
            batch_size:int=0, 
            folder:str="./models", 
            inference:bool=False, 
            **kwargs):

        super().__init__()
        self.data = {}
        self.data["dataname"] = name

        self.data["folder"] = folder
        self.data["epoch"] = None
        self.data["path"] = False
        self.data["model"] = self.__class__.__name__
        self.data["batch_size"] = batch_size
        self.mode = "default"
        self.writer = None
        self.mounted = False
        self.data["in_shape"] = in_shape
        self.inference = inference
        os.makedirs(self.data["folder"], exist_ok=True)

        self.criterions = {
            "default": nn.MSELoss()
        }
        self.optimisers = {
            "default": {
                "layer_set": "default",
                "class":torch.optim.Adam,
                "kwargs":{"lr":lr},
            }
        }
        self.schedulers = {
            "default": None
        }

        self.layers = {
            "default": {
            }
        }
        self.submodels = {
        }
        self.running_loss = {"total":0, "default":0}
        self.batch_loss = {"current_n":0, "current_list":[], "cumulative":0, "n_batches": 0}

        self.data["name"] = str(self)

        log("header", f"Model initialised: {self.data["name"]}")


    def __str__(self):
        if self.mounted:
            return f"{self.__class__.__name__}_{self.data['dataname']}"
        else:
            return f"{self.__class__.__name__}_{self.data['dataname']}"


    def __repr__(self):

        try: optim = self.optimisers[self.mode]
        except KeyError: optim = self.optimisers['default']
        try: crit = self.criterions[self.mode]
        except KeyError: crit = self.criterions['default']
        try: layers = self.layers[self.mode]
        except KeyError: layers = self.layers['default']
        try: loss = self.running_loss[self.mode]
        except KeyError: loss = self.running_loss['default']

        return f"<bi.CustomModel: {self.data['name']}\n - MODE: {self.mode}\n - class: {self.__class__.__name__}\n - optimiser: {optim.__class__.__name__}\n - criterion: {crit.__class__.__name__}\n - current epoch: {self.data['epoch']}\n - running loss: {loss}\n - layers: {[f'{k}({v.__class__.__name__})' for k, v in layers.items()]}\n>\n"


    def json(self):
        return json.dumps(self.data, indent=4)

    def to(self, device=None):
        if not self.mounted:
            self._mount_submodels()
        if device is None:
            device = DEVICE
        for k in self.submodels.keys():
            self.submodels[k] = self.submodels[k].to(device)
        return self


    def mount(self):
        self._mount_submodels
        return self

    def _mount_submodels(self):
        log(1, "Mounting submodels...")
        for k, layer_set in self.layers.items():
            self.submodels[k] = nn.Sequential(*[l.to(DEVICE) for l in layer_set.values()]).to(DEVICE)
        if not self.inference:
            for o, op_data in self.optimisers.items():
                self.optimisers[o] = self.optimisers[o]["class"](self.submodels[self.optimisers[o]["layer_set"]].parameters(), **self.optimisers[o]["kwargs"])
                if op_data.get("LRS", None) is not None:
                    self.schedulers[o] = op_data["LRS"](optimiser=self.optimisers[o])

        self.mounted = True

        self.reset_loss()
        self.data["name"] = str(self)

        if self.writer is None:
            self._create_writer()
        #self.writer.add_graph(self, torch.rand(self.data["in_shape"]))


        print(repr(self))
        self.to(DEVICE)


        return self.submodels.keys()


    def set_mode(self, mode:str):
        print(f"Model mode: {self.mode} -> {mode}")
        self.mode = mode
        return self.mode


    def __call__(self, x, **kwargs):
        return self.forward(x, **kwargs)


    def forward(self, x, submodel_name="mode"):
        if not self.mounted:
            self._mount_submodels()
        if submodel_name is None: return None
        if submodel_name == "mode": submodel_name = self.mode

        if submodel_name == "all":
            for submodel in self.submodels.values():
                x = submodel(x)
        else:
            x = self.submodels[submodel_name](x)
        return x


    def _create_writer(self):
        self.writer = SummaryWriter(log_dir=f"runs/{self.__class__.__name__}/{str(self)}_{datetime.datetime.now().strftime("%m-%d_%H-%M-%S")}")


    def reset_loss(self):
        self.running_loss["total"] = 0
        for c in self.running_loss.keys():
            self.running_loss[c] = 0
        self.batch_loss = {"current_n":0, "current_list":[], "cumulative":0, "n_batches": 0}


    def send_run(self, *args, **kwargs):
        if self.writer is not None:
            try:
                folder = platform.node()
                #print(self.writer.__dict__)
                run_name = os.path.basename(self.writer.log_dir)
                print(run_name)
                file = os.path.join(self.writer.log_dir, os.listdir(self.writer.log_dir)[-1])
                self.writer.flush()
                return send_tensorboard_run(*args, folder=folder, run=run_name, file=file, **kwargs)
            except Exception as e:
                log("warning", "Error uploading run:", e)


    def write_loss(self):
        av_losses = {}
        print("writing losses")
        for c, rl in self.running_loss.items():
            if c == "total":
                pass
            print("  ", c)
            total = self.running_loss["total"]
            if isinstance(rl, torch.Tensor):
                rl = rl.item()
            print(total, c, rl)

            if total == 0: av_loss = rl
            else: av_loss = rl / total
            if self.writer is not None:
                print("writing", av_loss)
                self.writer.add_scalar(f"loss/{c}", float(av_loss), self.data["epoch"])
            av_losses[c] = av_loss


        if self.data["batch_size"] != 0:
            if self.batch_loss["n_batches"] == 0: av_batch_loss = self.batch_loss["cumulative"]
            else: av_batch_loss = self.batch_loss["cumulative"] / self.batch_loss["n_batches"]
            if self.writer is not None:
                print("writing", av_batch_loss)
                self.writer.add_scalar(f"loss/batch", float(av_batch_loss), self.data["epoch"])

        return av_losses



    def set_epoch(self, epoch):
        self.data["epoch"] = epoch
        return self.data["epoch"]

    def leftover_batch(self):
        if self.data["batch_size"] != 0:
            if self.batch_loss["current_n"] != 0:
                log(2, f"Leftover backpropagation... ({self.batch_loss['current_n']})")
                self.loss(force_backpropagation=True)

    def add_epoch(self):

        self.leftover_batch()

        if self.writer is not None:
            for set_name, layers in self.layers.items():
                for layer_name, layer in layers.items():
                    if hasattr(layer, "weight"):
                        self.writer.add_histogram(f"{set_name}/weight/{layer_name}", layer.weight, self.data["epoch"])
                    if hasattr(layer, "bias"):
                        self.writer.add_histogram(f"{set_name}/bias/{layer_name}", layer.bias,  self.data["epoch"])



        av_losses = self.write_loss()

        if self.data["epoch"] not in (None, 0):
            self.step_schedulers(running_loss=av_losses)
        self.reset_loss()
        if self.data["epoch"] is None: self.data["epoch"] = 1
        else: self.data["epoch"] += 1
        return self.data["epoch"]


    def get_fname(self, add_epoch=False) -> str:
        if add_epoch and self.data["epoch"] is not None:
            fname = f"{self.data['name']}_E{self.data['epoch']}"
        else:
            fname = f"{self.data['name']}"
        return fname


    def save(self, path=None, add_epoch=False, temp=False):
        log(1, f"Saving model (TEMP={temp})")
        if path is None:
            path = os.path.join(self.data["folder"], self.get_fname(add_epoch=add_epoch))
        for name, submodel in self.submodels.items():
            if not path.endswith(".model.pt"):
                model_path = path + f".{name}.model.pt"
            else:
               model_path = model_path.replace(".model.pt", f".{name}.model.pt") 
            if temp:
                model_path = model_path.replace(".model.pt", ".temp.model.pt")
            if name == "default":
                self.data["path"] = model_path
            torch.save(submodel.state_dict(), model_path)
        return self.export(path=path, add_epoch=add_epoch, temp=temp)


    def export(self, path=None, add_epoch=False, temp=False):
        log(1, f"Exporting model (TEMP={temp})")

        if path is None:
            path = os.path.join(self.data["folder"], self.get_fname(add_epoch=add_epoch))
        if path.endswith(".model.pt"):
            data_path = path.replace(".model.pt", ".data.json")
        elif not path.endswith(".data.json"):
            data_path = path + ".data.json"
        else:
            data_path = path
        if temp and not path.endswith(".temp.data.json"):
            data_path = data_path.replace(".data.json", ".temp.data.json")
        json.dump(self.data, open(data_path, "w"), indent=4)
        return data_path


    def load(self, data_path=None, epoch=None, weights_only=False):
        log(1, f"Loading model (weights_only={weights_only})")

        if data_path is None:
            data_path = os.path.join(self.data["folder"], self.get_fname(add_epoch=epoch))+".data.json"
        if not os.path.exists(data_path):
            raise ModelNotFound(data_path)
        raw_data = json.load(open(data_path, "r"))
        self.data = self.data | raw_data

        base_path = data_path.replace(".data.json", "")
        is_tmp = base_path.endswith(".temp")
        if is_tmp:
            base_path = base_path.replace(".temp", "")

        print("BASE", base_path)

        for name, submodel in self.submodels.items():
            model_path = f"{base_path}.{name}"
            if is_tmp:
                model_path += ".temp"
            model_path += ".model.pt"
            submodel.load_state_dict(torch.load(model_path, weights_only=weights_only, map_location=torch.device("cpu")))
        return self


    def add_map(self, dataset):
        self.data["label_to_index"] = dataset.data["label_to_index"]
        self.data["index_to_label"] = dataset.data["index_to_label"]
        if dataset.data.get("lab_count", False):
            self.data["lab_count"] = dataset.data["lab_count"]

        return self.data["label_to_index"], dataset.data["index_to_label"]

    def add_histogram(self, name, data):
        if self.writer is None:
            self._create_writer()
        if not(isinstance(data, np.ndarray) or isinstance(data, torch.Tensor)):
            data = np.array(data)
        self.writer.add_histogram(name, data,  self.data["epoch"])

    def add_hparams(self, hparams={}, hmetrics={}):
        self.writer.add_hparams(hparams, hmetrics, run_name=".")

    def add_protein(self, name, coords, colours):
        self.writer.add_mesh(name, coords, colours, global_step=self.data["epoch"])

    def add_text(self, name, text):
        if self.writer is None:
            self._create_writer()
        self.writer.add_text(name, text)

    def test(self, dataset, re_load=False, temp=False):
        self.leftover_batch()

        log("header", "Model Validation")
        if re_load:
            log(1, f"reloading saved model: {self.data['path']}")

            self.load()
        dataset.test()
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
                item.to(DEVICE)
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

        self.writer.add_scalar(f"accuracy/total", (correct / total) * 100, self.data["epoch"])
        #self.writer.add_scalar(f"accuracy/dual/contactability", (confusion["contactability"]["right"] / total) * 100, self.data["epoch"])
        #self.writer.add_scalar(f"accuracy/dual/outer", (confusion["outer"]["right"] / total) * 100, self.data["epoch"])


        weighted_accuracy = sum([( v[k]/sum(v.values()) ) * ( 1-(truths.count(k) / len(truths)) ) for k, v in confusion.items()])/len(set(truths))
        self.writer.add_scalar(f"accuracy/weighted", weighted_accuracy * 100, self.data["epoch"])




        log(1, f"Results: EPOCH:{self.data['epoch']-1} correct={correct}, total={total}, accuracy={(correct / total) * 100:2.3f}%, weighted={weighted_accuracy * 100:2.3f}% ")

        if len(preds) > 0 and len(truths) > 0:
            try:
                from ..visualisation.plots import plot_confusion
                _, confusion_path = plot_confusion(preds, truths, title=f"{str(self)}", classes = label_to_index.keys())
                _, confusion_path = plot_confusion(preds, truths, title=f"{str(self)}.weighted", classes = label_to_index.keys())


                im = PIL.Image.open(confusion_path)
                im = torchvision.transforms.v2.functional.pil_to_tensor(im)
                self.writer.add_image(f"confusion", im, self.data["epoch"])
                del im
            except Exception as e:
                print(e)

        self.save(temp=temp)




    def _calculate_loss(self, output:torch.Tensor, item:Item, criterion_name:str="mode") -> torch.Tensor|float:

        if criterion_name == "mode": criterions = [self.mode]

        elif criterion_name == "all":
            criterions = [n for n in self.criterions.keys()]
        else:
            criterions = [criterion_name]


        losses = []

        for n, criterion in enumerate(criterions):
            if criterion not in self.criterions: criterions[n] = "default"; criterion = "default"
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

        if self.data["batch_size"] != 0:
            self.batch_loss["current_n"] += 1
            self.batch_loss["current_list"].append(loss)


        self.running_loss["total"] += 1
        #print(criterions)
        for l, c in zip(losses, criterions):
            self.running_loss[c] += l.item()

        return loss


    def _backpropagate(self, loss:torch.Tensor|None=None, zero_optims:str|None="mode", step="mode", force=False):
        if self.data["batch_size"] == 0:
            assert loss is not None
            loss.backward()
            self.step(step)
            self.zero_grad(zero_optims)
        else:
            #loss.backward(retain_graph=True)

            if self.batch_loss["current_n"] >= self.data["batch_size"] or force:
                batch_loss = torch.mean(torch.stack(self.batch_loss["current_list"]))
                print(f"Backpropagating ({batch_loss:5.4f})   ", end="\r")

                batch_loss.backward()

                self.batch_loss["n_batches"] += 1
                self.batch_loss["cumulative"] += batch_loss.item()

                self.batch_loss["current_n"] = 0
                self.batch_loss["current_list"] = []
                self.step(step)
                self.zero_grad(zero_optims)
        return loss





    def loss(self, output:torch.Tensor|None=None, item:Item|None=None, criterion_name:str="mode", backwards:bool=True, zero_optims:str|None="mode", step="mode", force_backpropagation=False) -> torch.Tensor|float:

        if (output is not None) and (item is not None):
            loss = self._calculate_loss(output, item, criterion_name=criterion_name)
        else: loss = None
        loss = self._backpropagate(loss=loss, zero_optims=zero_optims, step=step, force=force_backpropagation)

        return loss



    def step_schedulers(self, scheduler_name:str|None="mode", running_loss=None) -> bool:
        if scheduler_name is None: return False
        if scheduler_name == "mode": scheduler_name = self.mode
        if scheduler_name not in self.schedulers: scheduler_name = "default"

        if running_loss is None:
            running_loss = self.running_loss.get("default", 0.5)

        if scheduler_name == "all":
            for name, scheduler in self.schedulers.items():
                if scheduler is not None:
                    scheduler.step(running_loss=running_loss.get(name, 0.5))
        else:
            if self.schedulers[scheduler_name] is not None:
                lr = self.schedulers[scheduler_name].step(running_loss=running_loss.get(scheduler_name, 0.5))
                self.writer.add_scalar(f"learning_rate/{scheduler_name}", torch.Tensor(lr), self.data["epoch"])

        return True


    def step(self, optimizer_name:str|None="mode") -> bool:
        if optimizer_name is None: return False
        if optimizer_name == "mode": optimizer_name = self.mode
        if optimizer_name not in self.optimisers: optimizer_name = "default"


        if optimizer_name == "all":
            for optimizer in self.optimisers.values():
                optimizer.step()
        else:
            self.optimisers[optimizer_name].step()
        return True


    def zero_grad(self, optimizer_name:str|None="mode") -> bool:
        if optimizer_name is None: return False
        if optimizer_name == "mode": optimizer_name = self.mode
        if optimizer_name not in self.optimisers: optimizer_name = "default"

        if optimizer_name == "all":
            for optimizer in self.optimisers.values():
                optimizer.zero_grad()
        else:
            self.optimisers[optimizer_name].zero_grad()
        return True




