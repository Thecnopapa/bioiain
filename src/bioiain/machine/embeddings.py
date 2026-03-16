import os, json


from ..utilities.logging import log
from ..utilities.parallel import avail_cpus
from ..utilities.sequences import d3to1 ,d3toint

import torch


device = "cpu"


class MissingProgram(Exception):
    pass



class FoldseekError(Exception):
    pass



class Embedding(object):
    def __init__(self, *args, name=None, folder=None, **kwargs):
        assert name is not None
        self.name=name
        if folder is None:
            folder = "./embeddings"
        self.folder = os.path.join(folder, self.__class__.__name__+"s")
        self.subfolder = os.path.join(self.folder, self.name)
        self.path = os.path.join(self.subfolder, f"{self.name}.embedding.pt")
        os.makedirs(self.subfolder, exist_ok=True)
        self.length = 0
        self.iter_dim = 1


    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={self.length} at: {self.path}"


    def from_file(self, path, iter_dim=0):
        tensor = torch.load(path)
        self.name = path.split(".")[0]
        self.path = path
        self.folder = os.path.dirname(path)
        self.length = tensor.shape[iter_dim]
        self.iter_dim = iter_dim


    def get_tensor(self):
        tensor = torch.load(self.path)
        return tensor


    def generate_embedding(self, *args, **kwargs):
        raise NotImplementedError("Embedding: generate_embedding() must be overridden by subclass")



class PerResidueEmbedding(Embedding):
    def __init__(self, *args, entity=None, **kwargs):
        assert entity is not None
        super().__init__(self, *args, name=entity.name(), **kwargs)
        self.entity = entity
        self.subfolder = os.path.join(self.folder, self.name)

        self.sequence = self.entity.sequence()



class CVEmbedding(PerResidueEmbedding):
    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.entity = self.entity.fragment()
        self.cvectors = self.entity.cvectors()
        self.cvmatrix = self.entity.cvmatrix()


    def generate_embedding(self, *args, modulo_norm=2.3, **kwargs):

        e = []
        seq = ""
        for cv in self.cvectors:
            i = cv
            j = cv.closest
            i_j = cv.closest_vp

            len_i = i.d / modulo_norm
            len_j = j.d / modulo_norm
            len_i_j = i_j.d
            angle_i_j = i_j.a

            ri = d3toint[cv.resname]
            rn = d3to1[cv.resname]
            seq += rn

            e.append([len_i, len_j, angle_i_j, len_i_j, ri])

        e = torch.Tensor(e)
        torch.save(e, self.path)
        print(e, e.shape, len(seq))
        self.sequence = seq
        self.length = len(self.sequence)
        return self.path





class SaProtEmbedding(PerResidueEmbedding):
    def __init__(self, *args, foldseek_cmd="foldseek", with_foldseek=True, force=False, padding=1, keep_padding=False, use_temp=True, download_folder=".models", **kwargs):
        super().__init__(self, *args, **kwargs)
        self.folder = os.path.join(self.folder, "SaProt")
        self.subfolder = os.path.join(self.folder, self.name)
        self.fs_tokens = None
        self.foldseek_cmd = foldseek_cmd
        self.single_file = True
        self.length = len(self.sequence)
        self.keep_padding = keep_padding
        self.padding = padding
        self.download_folder = download_folder


        if keep_padding:
            self.length += padding*2



        self.iter_dim = 1
        self.generate_embedding(with_foldseek=with_foldseek, force=force, use_temp=use_temp)


    def generate_embedding(self, *args, with_foldseek=True, force=False, use_temp=False, **kwargs):
        #print("GENERATING_EMBEDDING")
        if with_foldseek:
            if self._get_foldseek(force=force, use_temp=use_temp) is None: return None
        return self._get_saprot(force=force)


    def _run_foldseek(self, out_path):
        #print("RUNNING FOLDSEEK")
        import subprocess
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        cmd = [self.foldseek_cmd, "structureto3didescriptor", "-v", "0", "--threads", f"{avail_cpus}", "--chain-name-mode", "0",
               self.entity.paths["self"], out_path]
        log("debug", "$", " ".join(cmd))
        subprocess.run(cmd)
        log("debug", ">", out_path)
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
                print("foldseek path: open", out_path)
                print(f.read())
                raise FoldseekError(f"No Foldseek data for: {self.entity}")
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


    def _get_foldseek(self, force=False, use_temp=True):
        #print("GETTING_FOLDSEEK")
        if use_temp:
            out_path = f"/tmp/bioiain/foldseek/{self.name}.foldseek.tsv"
        else:
            out_path = f"./foldseek/{self.name}.foldseek.tsv"

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
            #print("USING PRECALCULATED SAPROT at:",save_path)
            self.path = save_path
            return self

        if self.fs_tokens is None:
            tokenizer_name = "westlake-repl/SaProt_650M_PDB"
            model_name = "westlake-repl/SaProt_650M_PDB"
            in_tokens = [f"{s}#" for s in self.sequence]
        else:
            tokenizer_name = "westlake-repl/SaProt_650M_PDB"
            model_name = "westlake-repl/SaProt_650M_PDB"
            try:
                assert len(self.sequence) == len(self.fs_tokens)
            except AssertionError as e:
                print(self.sequence, len(self.sequence))
                print(self.fs_tokens, len(self.fs_tokens))
                raise FoldseekError("Foldseek length does not match sequence length")

            in_tokens = [f"{s.upper()}{fs.lower()}" for s, fs in zip(self.sequence, self.fs_tokens)]

        model_folder = self.download_folder


        #stop_logging()
        tokenizer_path = os.path.join(model_folder, f"tok_{tokenizer_name}")
        print("tok", tokenizer_path)
        if not os.path.exists(tokenizer_path):
            print("Downloading tokeniser...")
            tokenizer = AutoTokenizer.from_pretrained(tokenizer_name)
            os.makedirs(os.path.dirname(tokenizer_path), exist_ok=True)
            tokenizer.save_pretrained(tokenizer_path)
        tokenizer = AutoTokenizer.from_pretrained(tokenizer_path)

        model_path = os.path.join(model_folder, f"mod_{model_name}")
        print("mod", model_path)

        if not os.path.exists(model_path):
            print("Downloading model...")
            model = AutoModelForMaskedLM.from_pretrained(model_name)
            os.makedirs(os.path.dirname(model_path), exist_ok=True)
            model.save_pretrained(model_path)
        model = AutoModelForMaskedLM.from_pretrained(model_path)

        print("model eval")
        model.eval()
        model.to(device)
        #resume_logging()

        print("Preparing tokens...")
        inputs = tokenizer("".join(in_tokens), return_tensors="pt").to(device)
        inputs = {k: v.to(device) for k, v in inputs.items()}
        print("Tokens ready")
        #print(inputs)

        with torch.no_grad():
            print("Generating output...")
            outputs = model(**inputs, output_hidden_states=True)
            print("Output ready")

        last_hidden = outputs.hidden_states[-1]
        if not self.keep_padding:
            print(last_hidden.shape)
            print(self.keep_padding, self.padding)
            last_hidden = last_hidden[:,self.padding:-self.padding,:]
            print(last_hidden.shape)

        try:
            assert last_hidden.shape[1] == self.length
        except Exception as e:
            print(self)
            print(in_tokens)
            print(len(in_tokens))
            print(inputs)
            print(len(inputs))
            print(last_hidden)
            print(last_hidden.shape)
            print(self.sequence, self.length)
            print(last_hidden.shape[1], self.length)

            raise e
        torch.save(last_hidden, save_path)
        self.path = save_path
        #print("EMBEDDING SAVED AT:")
        #print(self.path)

        return self.path





















