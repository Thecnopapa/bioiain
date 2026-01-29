import os, json







#BASE CLASSES

class Embedding(object):

    def __init__(self, *args, name=None, folder=None, **kwargs):
        if folder is None:
            fodler = "./embeddings"
        self.folder =folder
        self.name=name
        self.path = None

    def from_file(self, path):
        self.name = path.split(".")[0]
        self.path = path
        self.folder = os.path.dirname(path)





class EmbeddingList(object):
    def __init__(self,*args,  name, folder="./embeddings", **kwargs):
        self.name = name
        self.folder = folder
        self.embeddings = {}


    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} N={len(self.embeddings)}>"

    def generate_embeddings(self, *args, **kwargs):
        raise NotImplementedError("EmbeddingList: generate_embeddings() must be overridden by subclass")


    def add(self, embedding:Embedding, key:str|int|None=None, label=None):
        if key is None:
            key = str(len(self.embeddings))

        self.embeddings[key] = {
                "key": key,
                "embedding_path": embedding.path,
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
            "embeddint_class": self[0]["embedding"].__class__.__name__,
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

class MonomerEmbedding(Embedding):
    pass



#CUSTOM CLASSES


class SaProtEmbeddings(PerResidueEmbeddings):

    def __init__(self, *args, **kwargs):
        super().__init__(self, *args, **kwargs)
        self.folder = os.path.join(self.folder, "SaProt")
        self.fs_tokens = None


    def generate_embeddings(self, *args, **kwargs):
        self._get_foldseek()
        self._run_saprot(self.sequence, self.fs_tokens)




    def _run_saprot(self, sequence, fs_tokens):
        from transformers import AutoTokenizer, AutoModelForMaskedLM
        import torch

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

    def _get_foldseek(self):

    # def run_foldseek(filename, data_folder, raw_folder, label_folder):
    #     os.makedirs(raw_folder, exist_ok=True)
    #
    #     if os.path.exists(f"{raw_folder}/{filename}.foldseek.json") and not config["general"]["force"]:
    #         bi.log(3, "Foldseek already calculated")
    #         return True
    #     label_dict = json.load(open(f"{label_folder}/{filename}.labels.json"))
    #     cmd = [config["general"]["foldseek"],
    #            "structureto3didescriptor", "-v", "0", "--threads", "4",
    #            "--chain-name-mode", "0", f"{data_folder}/{filename}.cif",
    #            f"{raw_folder}/{filename}.foldseek.csv"
    #            ]
    #     bi.log(4, " ".join(cmd))
    #     subprocess.run(cmd)
    #     done_chains = []
    #     foldseek_dict = {k: None for k in label_dict.keys()}
    #     with open(f"{raw_folder}/{filename}.foldseek.csv", "r", encoding="utf-8") as f:
    #
    #         for line, ch in zip(f, foldseek_dict.keys()):
    #
    #             if ch in done_chains:
    #                 continue
    #             done_chains.append(ch)
    #             # print(ch)
    #
    #             rns, tks = line.split("\t")[1:3]
    #             resns = [r for r in rns]
    #             toks = [t for t in tks]
    #
    #             bi.log(3, "foldseek out:", ch, len(resns), len(toks))
    #
    #             if not len(resns) == len(toks):
    #                 print(resns)
    #                 print(toks)
    #                 print(len(resns), len(toks))
    #                 foldseek_dict.pop(ch)
    #                 bi.log("warning", "foldseek tokens and dssp_dict do not match:", ch)
    #                 exit()
    #             foldseek_dict[ch] = {}
    #             for n, (r, t) in enumerate(zip(resns, toks)):
    #                 if toks[n] == " ":
    #                     toks = "-"
    #                 foldseek_dict[ch][n] = {"fs": t, "resn": r}
    #     json.dump(foldseek_dict, open(f"{raw_folder}/{filename}.foldseek.json", "w"), indent=4)
    #     return True
            pass















