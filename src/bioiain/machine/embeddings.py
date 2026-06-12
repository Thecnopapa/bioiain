import os, json

from ..utilities import *
from ..utilities.parallel import avail_cpus
from ..utilities.sequences import d3
from ..utilities.exceptions import *

import torch
from torch import Tensor

import numpy as np
from typing_extensions import Self


device = "cpu"




class Embedding(object):
    param_names = []
    def __init__(self, name=None, folder=None, subfolder=None,group_by_class=True, **kwargs):
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

    def exists(self):
        return os.path.exists(self.path())

    def path(self, force=False) -> str:
        if self._path is not None and not force:
            return self._path
        path = self.folder
        if self.group_by_class:
            path = os.path.join(path, self.__class__.__name__)
        os.makedirs(path, exist_ok=True)
        if self.subfolder is not None:
            path = os.path.join(path, self.subfolder)
        path = os.path.join(path, self.name+".pt")
        self._path = path
        return path

    def tensor(self, force=False, generate=True) -> Tensor|None:
        if self._tensor is not None and not force:
            return self._tensor

        if self.exists():
            tensor = torch.load(self.path())
        elif generate:
            tensor = self.generate()
        else:
            return None
        self._tensor = tensor
        return tensor

    def save(self, **kwargs) -> Self:
        tensor = self.tensor(**kwargs)
        if tensor is not None:
            os.makedirs(os.path.dirname(self.path()), exist_ok=True)
            torch.save(tensor, self.path())
            return self
        else:
            raise EmptyTensor()

    def length(self) -> int:
        return len(self.tensor())

    def __len__(self):
        return self.length()

    @classmethod
    def from_file(cls, path, **kwargs):
        self = cls(**kwargs)
        self._path = path
        return self

    @classmethod
    def from_tensor(cls, tensor, **kwargs):
        self = cls(**kwargs)
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


class ProteinEmbedding(Embedding):
    residue_embedding_class = None
    def __init__(self, entity=None, residue_embedding_class=None, **kwargs):
        super().__init__(**kwargs)
        self.entity = entity
        self.sequence = None
        if residue_embedding_class is not None:
            self.residue_embedding_class = residue_embedding_class
            self.param_names = self.residue_embedding_class.param_names
        self.missing_indexes = []
        if self.entity is not None and not self.entity.has_flag("no_atoms", True):
            if self.name == self.__class__.__name__:
                self.name = self.entity.name()
                self.entity_path = self.entity.path()
            self.sequence = self.entity.sequence()

    def _generate(self, *args, **kwargs) -> list:
        assert self.residue_embedding_class is not None
        e = []
        for n, res in enumerate(self.entity.residues()):
            try:
                e.append(self.residue_embedding_class(*args, residue=res,  **kwargs).tensor())
            except NoEmbeddingForThisResidue:
                self.missing_indexes.append(n)
        return e























