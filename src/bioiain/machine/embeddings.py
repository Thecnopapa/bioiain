import os, json


from ..utilities.logging import log



device = "cpu"


#BASE CLASSES

class Embedding(object):

    def __init__(self, *args, name=None, folder=None, **kwargs):
        if folder is None:
            folder = "./embeddings"
        self.folder =folder
        os.makedirs(self.folder, exist_ok=True)
        self.name=name
        self.path = None
        self.has_multiple = None
        self.range=None

    def from_file(self, path, has_multiple=False, target_range=(None,None)):
        self.name = path.split(".")[0]
        self.path = path
        self.folder = os.path.dirname(path)
        if has_multiple:
            self.has_multiple = True
            self.range = target_range





class EmbeddingList(object):
    def __init__(self,*args,  name, folder="./embeddings", single_file=False, **kwargs):
        self.name = name
        self.folder = folder
        os.makedirs(self.folder, exist_ok=True)
        self.embeddings = {}
        self.single_file = single_file


    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={len(self.embeddings)}>"

    def generate_embeddings(self, *args, **kwargs):
        raise NotImplementedError("EmbeddingList: generate_embeddings() must be overridden by subclass")


    def add(self, embedding:Embedding, key:str|int|None=None, label=None):
        if key is None:
            key = str(len(self.embeddings))

        if embedding.has_multiple:
            self.embeddings[key] = {
            "range": embedding.range,
            "embedding_path": embedding.path,
            "label_path": label,
            }

        else:
            self.embeddings[key] = {
                "key": key,
                "embedding_path":embedding.path,
                "label": label,
            }
        return self[key]

    def __getitem__(self, key):
        return self.embeddings[key]

    def add_label(self, key, label):
        self.embeddings[key]["label"] = label
        return self[key]

    def export(self, folder=None):
        if folder is None:
            assert self.folder is not None
            folder = self.folder
        data = {
            "name": self.name,
            "embedding_class": self[0]["embedding"].__class__.__name__,
            "list_class": self.__class__.__name__,
            "embeddings": self.embeddings,
        }
        fname = f"{self.name}.{data['list_class']}.embeddings.json"
        path = os.path.join(folder, fname)
        json.dump(data, open(path, "w"))
        return path



class ResidueEmbedding(Embedding):
    pass



class PerResidueEmbeddings(EmbeddingList):
    def __init__(self, *args, entity, **kwargs):
        super().__init__(self, *args, name=entity.get_name(), **kwargs)
        self.sequence = None
        self.entity = entity
        self._get_sequence()

    def _get_sequence(self):
        self.sequence = self.entity.get_sequence(True)
        return self.sequence

class MonomerEmbedding(Embedding):
    pass



#CUSTOM CLASSES


class SaProtEmbeddings(PerResidueEmbeddings):

    def __init__(self, *args, foldseek_cmd="foldseek", **kwargs):
        super().__init__(self, *args, **kwargs)
        self.folder = os.path.join(self.folder, "SaProt")
        os.makedirs(self.folder, exist_ok=True)
        self.fs_tokens = None
        self.foldseek_cmd = foldseek_cmd
        self.single_file = True


    def generate_embeddings(self, *args, with_foldseek=True, **kwargs):
        if with_foldseek:
            self._get_foldseek()
        self._run_saprot()


    def _get_foldseek(self, force=False):
        import subprocess
        out_path = f"/tmp/bioiain/foldseek/{self.name}.foldseek.tsv"
        if not os.path.exists(out_path) or force:
            os.makedirs(os.path.dirname(out_path), exist_ok=True)
            cmd = [self.foldseek_cmd, "structureto3didescriptor", "-v", "0", "--threads", "4", "--chain-name-mode", "0", self.entity.paths["self"], out_path]
            log("debug", " ".join(cmd))
            subprocess.run(cmd)

        with open(out_path, "r", encoding="utf-8") as f:
            raw = f.read().split("\t")
            fname, seq, tokens = raw[:3]
            assert seq.strip() == self.sequence
            self.fs_tokens = tokens.strip()
            return self.fs_tokens
        return None


    def _run_saprot(self):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch


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
        #print(last_hidden.shape, len(self.sequence))
        save_path = os.path.join(self.folder, f"{self.name}.embedding.pt")
        torch.save(last_hidden, save_path)

        embedding = Embedding(name=self.name)
        embedding.from_file(save_path, has_multiple=True)
        self.add(embedding)
        return embedding






        # def run_saprot(name, mode, foldseek_path, label_path, save_folder):
    #     bi.log(3, "Running SaProt, mode:", mode)
    #     label_dict = json.load(open(f"{label_path}/{name}.labels.json"))
    #     # print(label_dict.keys())
    #     fs_keys = label_dict.keys()
    #     if mode == "full":
    #         assert foldseek_path is not None
    #         foldsek_dict = json.load(open(f"{foldseek_path}/{name}.foldseek.json"))
    #         # print(foldsek_dict.keys())
    #         assert label_dict.keys() == foldsek_dict.keys()
    #         fs_keys = foldsek_dict.keys()
    #
    #     seqs = {}
    #     # print(label_dict.keys(), foldsek_dict.keys())
    #
    #     for ch, fch in zip(label_dict.keys(), fs_keys):
    #         bi.log(4, "Merging foldseek_dict:", ch, fch)
    #         if mode == "full":
    #             if foldsek_dict[ch] is None:
    #                 bi.log("warning", f"chain {ch} has no foldseek data")
    #                 continue
    #             # print(foldsek_dict[ch])
    #
    #             if len(label_dict[ch]) != len(foldsek_dict[ch]):
    #                 bi.log("warning", "label and foldseek_dict do not match:", ch, len(label_dict[ch]),
    #                        len(foldsek_dict[ch]))
    #                 continue
    #             try:
    #                 seqs[ch] = [f"{l['resn'].upper()}{f['fs'].lower()}" for l, f in
    #                             zip(label_dict[ch].values(), foldsek_dict[ch].values())]
    #             except:
    #                 bi.log("warning", "unknown atom in chain:", ch)
    #                 [bi.log("warning", f"{r['res']} -> {r['resn']} / {r['resn3']}") for r in label_dict[ch].values()
    #                  if None in [r["res"], r["resn"], r["resn3"]]]
    #         elif mode == "seq":
    #             seqs[ch] = [f"{l['resn']}#" for l in label_dict[ch].values()]
    #         else:
    #             bi.log("error", "Unknown SaProt mode:", mode)
    #     # print("FOLDSEEK", seqs.keys())
    #
    #     for ch in seqs.keys():
    #         # Load model directly
    #         bi.log(4, "Generating embeddings:", ch)
    #         device = config["general"]["device"]
    #         if mode == "full":
    #             tokenizer = AutoTokenizer.from_pretrained("westlake-repl/SaProt_35M_AF2")
    #             model = AutoModelForMaskedLM.from_pretrained("westlake-repl/SaProt_35M_AF2")
    #         elif mode == "seq":
    #             tokenizer = AutoTokenizer.from_pretrained("westlake-repl/SaProt_35M_AF2_seqOnly")
    #             model = AutoModelForMaskedLM.from_pretrained("westlake-repl/SaProt_35M_AF2_seqOnly")
    #         else:
    #             bi.log("error", "Unknown SaProt mode:", mode)
    #
    #         model.eval()
    #         model.to(device)
    #
    #         seq = "".join(seqs[ch])
    #
    #         inputs = tokenizer(seq, return_tensors="pt").to(device)
    #         inputs = {k: v.to(device) for k, v in inputs.items()}
    #         # print(inputs)
    #
    #         with torch.no_grad():
    #             outputs = model(**inputs, output_hidden_states=True)
    #         # print(outputs)
    #
    #         # outputs.hidden_states is a tuple of all layers, including embeddings
    #         # Shape of each layer: [batch_size, sequence_length, hidden_dim]
    #         all_hidden_states = outputs.hidden_states
    #
    #         # Last layer hidden states
    #         last_hidden = all_hidden_states[-1]  # [1, seq_len, hidden_dim]
    #         # print(last_hidden.shape)  # ['<cls>', 'M#', 'E#', 'V#', 'Q#', '<eos>']
    #         # print(last_hidden)
    #
    #         os.makedirs(save_folder, exist_ok=True)
    #         torch.save(last_hidden, f"{save_folder}/{name}_{ch}.pt")
    #
    #     return True
        pass




















