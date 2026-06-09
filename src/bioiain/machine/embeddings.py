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
    def __init__(self, name=None, folder=None, subfolder=None, **kwargs):
        self.name = name
        self.folder = folder
        if self.folder is None:
            self.folder = os.path.join(SUBDIR_NAME, "embeddings")

        self.subfolder = subfolder
        self.length = None

        self._path = None
        self._tensor = None

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={self.length} at: {self.path()}>"

    def path(self, force=False) -> str:
        if self._path is not None and not force:
            return self._path
        path = self.folder
        if self.subfolder is not None:
            path = os.path.join(path, self.subfolder)
        path = os.path.join(path, self.name+".pt")
        self._path = path
        return path

    def tensor(self, force=False, generate=False) -> Tensor|None:
        if self._tensor is not None and not force:
            return self._tensor

        if os.path.exists(self.path()):
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
            torch.save(tensor, self.path())
            return self
        else:
            raise EmptyTensor()


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

    def _generate(self, *args, **kwargs) -> Tensor|list|np.array:
        raise NotImplementedError("Embedding: _generate() must be overridden by subclass")

    def generate(self, *args, append=False, append_dim=0, **kwargs) -> Tensor:
        t = self._generate(*args, **kwargs)
        if not isinstance(t, Tensor):
            t = Tensor(t)
        if self.tensor() is not None and append:
            t = torch.cat((self.tensor(), t), dim=append_dim)
        self._tensor = t
        return self._tensor


class ResidueEmbedding(Embedding):
    def __init__(self, residue=None, **kwargs):
        super().__init__(**kwargs)
        self.residue = residue
        if self.residue is not None and self.name is None:
            self.name = self.residue.name()
        self.param_names = []

    def get_param_name(self, pos, when_missing=None):
        try:
            return self.param_names[pos]
        except IndexError:
            return when_missing


class ProteinEmbedding(Embedding):
    def __init__(self, entity=None, **kwargs):
        super().__init__(**kwargs)
        self.entity = entity
        self.residue_embedding_class = None
        self.residue_embeddings = []
        if self.entity is not None and self.name is None:
            self.name = self.entity.name()






















