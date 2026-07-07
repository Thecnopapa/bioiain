import os, json

from ..utilities import *
from ..utilities.sequences import d3
from ..utilities.exceptions import *
from ..utilities.files import relative_path

import torch
from torch import Tensor

import numpy as np
from typing_extensions import Self


device = "cpu"




class Embedding(object):
    param_names = []
    def __init__(self, name=None, folder=None, subfolder=None,group_by_class=True, dry=False, **kwargs):
        self.dry = dry
        if name is not None:
            self.name = name
        else:
            self.name = self.__class__.__name__
        self.folder = folder
        if self.folder is None:
            self.folder = os.path.join(SUBDIR_NAME, "embeddings")

        self.group_by_class = group_by_class
        self.subfolder = subfolder
        self._path = None
        self._tensor = None

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={self.length()} at: {self.path()}>"


    def dict(self, extra:dict={}):
        return {
            "name": self.name,
            "folder": relative_path(self.folder),
            "subfolder": relative_path(self.subfolder),
            "group_by_class": self.group_by_class,
            "embedding_path": relative_path(self.path()),
            "length": len(self),
            "iter_dim": getattr(self, "iter_dim", 0),
            "shape": self.tensor().shape,
            "param_names": getattr(self, "param_names", None),
        } | extra


    def exists(self, check_json=True):
        if check_json:
            return os.path.exists(self.json()) and os.path.exists(self.path())
        return os.path.exists(self.path())

    def json(self):
        return self.path().replace(".pt", ".json")

    def reload(self):
        return self.__class__.from_json(self.json())

    def path(self, force=False) -> str:
        if self._path is not None and not force:
            return self._path
        path = self.folder
        if self.group_by_class:
            path = os.path.join(path, self.__class__.__name__)
            if getattr(self, "residue_embedding_name", None) is not None:
                path = os.path.join(path, self.residue_embedding_name)
        os.makedirs(path, exist_ok=True)
        if self.subfolder is not None:
            path = os.path.join(path, self.subfolder)
        path = os.path.join(path, self.name+".pt")
        self._path = path
        return path

    def tensor(self, force=False, generate=True) -> Tensor|None:
        if self._tensor is not None and not force:
            #print("Tensor cached")
            return self._tensor

        if self.exists():
            tensor = torch.load(self.path())
        elif generate:
            tensor = self.generate()
        else:
            raise NoTensorAvailable()
        self._tensor = tensor
        return tensor

    def save(self, **kwargs) -> Self:
        tensor = self.tensor(**kwargs)
        if tensor is not None:
            os.makedirs(os.path.dirname(self.path()), exist_ok=True)
            torch.save(tensor, self.path())
            json.dump(self.dict(), open(self.json(), "w"), indent=4)
            return self
        else:
            raise EmptyTensor()

    def length(self) -> int:
        if self.dry:
            return 0
        iter_dim = getattr(self, "iter_dim", 0)
        if iter_dim is not None:
            return self.tensor().shape[iter_dim]
        return 1

    def __len__(self):
        return self.length()

    @classmethod
    def from_file(cls, path, force_as_tensor=False, force_as_json=False, **kwargs):
        if path.endswith(".pt") or force_as_tensor:
            self = cls.from_tensor(path, **kwargs)
            self._path = path
        elif path.endswith(".json") or force_as_json:
            self = cls.from_json(path, **kwargs)
        else:
            raise UnknownEmbeddingFormat()
        return self


    @classmethod
    def from_json(cls, path, **kwargs):
        self = cls()
        data = json.load(open(path))
        for k, v in data.items():
            if k == "length":
                continue
            setattr(self, k, v)
        return self

    @classmethod
    def from_tensor(cls, tensor, **kwargs):
        self = cls(**kwargs)
        if type(tensor) is torch.Tensor:
            pass
        elif type(tensor) is str:
            tensor = torch.load(self.path())
        elif type(tensor) in (list, tuple, np.ndarray):
            tensor = torch.tensor(np.array(tensor))
        self._tensor = tensor
        return self

    def append(self, t, append_dim=0):
        if self.tensor() is not None:
            t = torch.cat((self.tensor(), t), dim=append_dim)
        self._tensor = t
        return self._tensor

    def _generate(self, *args, **kwargs) -> list:
        raise NotImplementedError("Embedding: _generate() must be overridden by subclass")

    def generate(self, *args, append=False, append_dim=0, **kwargs) -> Tensor:
        t = self._generate(*args, **kwargs)
        if not isinstance(t, Tensor):
            t = np.array(t)
            #print(t.shape)
            t = Tensor(t)
        if append:
            self.append(t, append_dim)
        else:
            self._tensor = t
        return self._tensor


class ResidueEmbedding(Embedding):
    def __init__(self, residue=None, **kwargs):
        super().__init__(**kwargs)
        self.residue = residue
        if self.residue is not None and self.name == self.__class__.__name__:
            self.name = self.residue.name()
        self.param_names = []

    def get_param_name(self, pos, when_missing=None):
        try:
            return self.param_names[pos]
        except IndexError:
            return when_missing

    def save(self):
        raise NotAGoodIdea()

    def dict(self, extra:dict={}):
        return super().dict({"residue": self.residue, "entity":str(self.param_names), "entity_path":self.entity_path}|extra)


class ProteinEmbedding(Embedding):
    residue_embedding_class = None
    residue_embedding_name = None
    def __init__(self, entity=None, residue_embedding_class=None, **kwargs):
        super().__init__(**kwargs)
        self.entity = entity
        self.sequence = None
        self.chains = "*"
        if residue_embedding_class is not None:
            self.residue_embedding_class = residue_embedding_class
            self.param_names = self.residue_embedding_class.param_names
            self.residue_embedding_name = self.residue_embedding_class.__name__

        self.missing_indexes = []


        if self.entity is not None :
            if self.name == self.__class__.__name__:
                self.name = self.entity.name()
                if not self.entity.has_flag("no_atoms", True):
                    self.entity_path = self.entity.path()

    def dict(self, extra={}):
        return super().dict({"sequence": self.sequence,
                             "entity":str(self.entity),
                             "entity_path":self.entity_path,
                             "chains":self.chains,
                             "missing_indexes":self.missing_indexes,
                             "residue_embedding_name":self.residue_embedding_class.__name__,
                             }|extra)

    def _generate(self, *args, **kwargs) -> list:
        assert self.residue_embedding_class is not None
        e = []
        seq = ""
        for n, res in enumerate(self.entity.residues()):
            try:
                e.append(self.residue_embedding_class(*args, residue=res,  **kwargs).tensor(force=True))
                seq += d3(res.resname)[0]
            except NoEmbeddingForThisResidue:
                seq += "-"
                self.missing_indexes.append(n)
        self.sequence = seq
        return e
