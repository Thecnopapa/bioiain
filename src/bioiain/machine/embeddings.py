import os, json


from ..utilities.logging import log

import torch
from torch.utils.data import Dataset

device = "cpu"



class Embedding(object):

    def __init__(self, *args, name, folder=None, subfolder=None **kwargs):
        self.name=name
        if folder is None:
            folder = "./embeddings"
        self.folder =folder
        os.makedirs(self.folder, exist_ok=True)
        if subfolder is None:
            subfolder = self.name
        self.path = None
        self.length = 1

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={len(self.embeddings


    def from_file(self, path, has_multiple:bool=False, length:int=1):
        self.name = path.split(".")[0]
        self.path = path
        self.folder = os.path.dirname(path)
        self.length = length

    def get_tensor(self):
        tensor = torch.load(self.path)
        return tensor

    def generate_embedding(self, *args, **kwargs):
        raise NotImplementedError("Embedding: generate_embedding() must be overridden by subclass")



class PerResidueEmbedding(Embedding):
    def __init(self, *args, entity, **kwargs):
        super().__init__(self, *args, name=entity.get_name(), **kwargs)
        self.entity = entity
        self.sequence = self._get_sequence()
        self.subfolder = os.path.join(self.folder, self.name)
        os.makedirs(self.subfolder)

    def _get_sequence(self):
        self.sequence = self.entity.get_sequence()
        return self.sequence


class SaProtEmbedding(PerResidueEmbedding):

    def __init__(self, *args, foldseek_cmd="foldseek", **kwargs):
        super().__init__(self, *args, **kwargs)
        self.folder = os.path.join(self.folder, "SaProt")
        os.makedirs(self.folder, exist_ok=True)
        self.fs_tokens = None
        self.foldseek_cmd = foldseek_cmd
        self.single_file = True
        self.length = len(self.sequence)+2

    def generate_embedding(self, *args, with_foldseek=True, force=False, **kwargs):
        if with_foldseek:
            if self._get_foldseek(force=force) is None: return None
        return self._get_saprot(force=force)

    def _run_foldseek(self, out_path):
        import subprocess
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        cmd = [self.foldseek_cmd, "structureto3didescriptor", "-v", "0", "--threads", "4", "--chain-name-mode", "0",
               self.entity.paths["self"], out_path]
        log("debug", "$", " ".join(cmd))
        subprocess.run(cmd)

    def _read_foldseek(self, out_path):
        with open(out_path, "r", encoding="utf-8") as f:
            raw = f.read().split("\t")
            try:
                fname, seq, tokens = raw[:3]
            except:
                print(out_path, ":")
                print(f.read())
                return None
            try:
                seq.strip() == self.sequence
            except AssertionError:
                print(self.name)
                print("seq:", seq.strip())
                print("self.sequence:", self.sequence)
                raise
            self.fs_tokens = tokens.strip()
            return self.fs_tokens

    def _get_foldseek(self, force=False):
        out_path = f"/tmp/bioiain/foldseek/{self.name}.foldseek.tsv"
        if not os.path.exists(out_path) or force:
            self._run_foldseek(out_path)
            return self._read_foldseek(out_path)
        else:
            return self._read_foldseek(out_path)


    def _get_saprot(self, force=False):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch

        save_path = os.path.join(self.subfolder, f"{self.name}.embedding.pt")
        embedding = Embedding(name=self.name)

        if os.path.exists(save_path) and not force:
            return embedding.from_file(save_path, has_multiple=True)

        if self.fs_tokens is None:
            tokenizer_name = "westlake-repl/SaProt_35M_AF2_seqOnly"
            model_name = "westlake-repl/SaProt_35M_AF2_seqOnly"
            in_tokens = [f"{s}#" for s in self.sequence]
        else:
            tokenizer_name = "westlake-repl/SaProt_35M_AF2"
            model_name = "westlake-repl/SaProt_35M_AF2"
            in_tokens = [f"{s.upper()}{fs.lower()}" for s, fs in zip(self.sequence, self.fs_tokens)]

        #print(in_tokens)
        seq = "".join(in_tokens)
        #print(seq)

        tokenizer_path = f"/tmp/bioiain/models/tok_{tokenizer_name}"
        if not os.path.exists(tokenizer_path):
            tokenizer = AutoTokenizer.from_pretrained(tokenizer_name)
            os.makedirs(os.path.dirname(tokenizer_path), exist_ok=True)
            tokenizer.save_pretrained(tokenizer_path)
        tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)

        model_path = f"/tmp/bioiain/models/mod_{model_name}"
        if not os.path.exists(model_path):
            model = AutoModelForMaskedLM.from_pretrained(model_name)
            os.makedirs(os.path.dirname(model_path))
            model.save_pretrained(model_path, exist_ok=True)
        model = AutoModelForMaskedLM.from_pretrained(model_path)


        model.eval()
        model.to(device)


        inputs = tokenizer("".join(in_tokens), return_tensors="pt").to(device)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        #print(inputs)

        with torch.no_grad():
            outputs = model(**inputs, output_hidden_states=True)

        last_hidden = outputs.hidden_states[-1]
        assert last_hidden.shape[1] == self.length

        torch.save(last_hidden, save_path)
        self.path = save_path

        return self.save_path


















