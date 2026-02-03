import os, json


from ..utilities.logging import log

import torch
from torch.utils.data import Dataset

device = "cpu"

class MissingProgram(Exception):
    pass

class Embedding(object):

    def __init__(self, *args, name=None, folder=None, **kwargs):
        assert name is not None
        self.name=name
        if folder is None:
            folder = "./embeddings"
        self.folder =folder
        os.makedirs(self.folder, exist_ok=True)
        self.path = None
        self.length = 0
        self.iter_dim = 1

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={self.length} at: {self.path}"


    def from_file(self, path, iter_dim=0):
        tensor = torch.load(path)
        self.name = path.split(".")[0]
        self.path = path
        self.folder = os.path.dirname(path)
        self.length = tensor.shape[length_dim]
        self.iter_dim = iter_dim

    def get_tensor(self):
        tensor = torch.load(self.path)
        return tensor

    def generate_embedding(self, *args, **kwargs):
        raise NotImplementedError("Embedding: generate_embedding() must be overridden by subclass")



class PerResidueEmbedding(Embedding):
    def __init__(self, *args, entity=None, **kwargs):
        assert entity is not None
        super().__init__(self, *args, name=entity.get_name(), **kwargs)
        self.entity = entity
        self.sequence = self._get_sequence()
        self.subfolder = os.path.join(self.folder, self.name)

    def _get_sequence(self):
        self.sequence = self.entity.get_sequence()
        return self.sequence


class SaProtEmbedding(PerResidueEmbedding):

    def __init__(self, *args, foldseek_cmd="foldseek", with_foldseek=True, force=False, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.folder = os.path.join(self.folder, "SaProt")
        self.subfolder = os.path.join(self.folder, self.name)
        self.fs_tokens = None
        self.foldseek_cmd = foldseek_cmd
        self.single_file = True
        self.length = len(self.sequence)+2
        self.iter_dim = 1
        self.generate_embedding(with_foldseek=with_foldseek, force=force)

    def generate_embedding(self, *args, with_foldseek=True, force=False, **kwargs):
        #print("GENERATING_EMBEDDING")
        if with_foldseek:
            if self._get_foldseek(force=force) is None: return None
        return self._get_saprot(force=force)

    def _run_foldseek(self, out_path):
        #print("RUNNING FOLDSEEK")
        import subprocess
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        cmd = [self.foldseek_cmd, "structureto3didescriptor", "-v", "0", "--threads", "4", "--chain-name-mode", "0",
               self.entity.paths["self"], out_path]
        log("debug", "$", " ".join(cmd))
        subprocess.run(cmd)
        if not os.path.exists(out_path):
            raise MissingProgram("Foldseek not installed or not working")

    def _read_foldseek(self, out_path, try_again=False):
        #print("READING FOLDSEEK")
        with open(out_path, "r", encoding="utf-8") as f:
            raw = f.read().split("\t")
            try:
                fname, seq, tokens = raw[:3]
            except:
                if try_again:
                    self._run_foldseek(out_path)
                print(out_path, ":")
                print(f.read())
                raise Exception("No Foldseek data")
            try:
                seq.strip() == self.sequence
            except AssertionError:
                print(self.name)
                print("seq:", seq.strip())
                print("self.sequence:", self.sequence)
                raise
            self.fs_tokens = tokens.strip()
            #print(self.fs_tokens)
            return self.fs_tokens

    def _get_foldseek(self, force=False):
        #print("GETTING_FOLDSEEK")
        out_path = f"/tmp/bioiain/foldseek/{self.name}.foldseek.tsv"
        #print(not os.path.exists(out_path), force)
        if (not os.path.exists(out_path)) or force:
            self._run_foldseek(out_path)
            return self._read_foldseek(out_path)
        else:
            return self._read_foldseek(out_path, try_again=True)

    def _get_saprot(self, force=False):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch
        #print("GETTING_SAPROT")
        os.makedirs(self.subfolder, exist_ok=True)
        save_path = os.path.join(self.subfolder, f"{self.name}.embedding.pt")


        if os.path.exists(save_path) and not force:
            print("USING PRECALCULATED SAPROT at:",save_path)
            self.path = save_path
            return self

        if self.fs_tokens is None:
            tokenizer_name = "westlake-repl/SaProt_650M_PDB"
            model_name = "westlake-repl/SaProt_650M_PDB"
            in_tokens = [f"{s}#" for s in self.sequence]
        else:
            tokenizer_name = "westlake-repl/SaProt_650M_PDB"
            model_name = "westlake-repl/SaProt_650M_PDB"
            in_tokens = [f"{s.upper()}{fs.lower()}" for s, fs in zip(self.sequence, self.fs_tokens)]



        tokenizer_path = f"/tmp/bioiain/models/tok_{tokenizer_name}"
        if not os.path.exists(tokenizer_path):
            tokenizer = AutoTokenizer.from_pretrained(tokenizer_name)
            os.makedirs(os.path.dirname(tokenizer_path), exist_ok=True)
            tokenizer.save_pretrained(tokenizer_path)
        tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)

        model_path = f"/tmp/bioiain/models/mod_{model_name}"
        if not os.path.exists(model_path):
            model = AutoModelForMaskedLM.from_pretrained(model_name)
            os.makedirs(os.path.dirname(model_path), exist_ok=True)
            model.save_pretrained(model_path)
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
        print("EMBEDDING SAVED AT:")
        print(self.path)

        return self.path


















