import os, sys, json, subprocess

sys.path.append('..')

from src.bioiain.base import BIEntity, BIChain
from src.bioiain.aleph import FragmentedStructure
import polars as pl
import numpy as np
from src.bioiain.utilities.exceptions import *
from src.bioiain import log
from src.bioiain.utilities import *
from src.bioiain.utilities.maths import *
from src.bioiain.utilities.kdtree import KDT
from src.bioiain.utilities.files import StructureDataset
from src.bioiain.utilities.sequences import FASTA










class FoldseekDB(object):
    def __init__(self,name, list_dir_or_dataset, folder=None, foldseek_command=None, force=False, verbose=2, saprot_model="SaProt_650M_PDB", dry=False):
        if foldseek_command is None:
            foldseek_command = os.environ.get("FOLDSEEK_PATH", "foldseek")
        self.foldseek_command = foldseek_command
        self.saprot_model=saprot_model

        self.verbose = str(verbose)
        if folder is None:
            folder = os.path.join(SUBDIR_NAME, "foldseek")
        self.folder = os.path.join(folder, name)
        self.db_path = os.path.join(self.folder, name)
        os.makedirs(self.folder, exist_ok=True)
        self.name = name
        self.list_dir_or_dataset = list_dir_or_dataset

        self.token_fasta_path = None
        self.sequence_fasta_path = None
        self.force = force

        if not dry:
            self.run()


        


    def run(self, force=None):
        if force is None:
            force = self.force
        self.create_db(force=force)
        self.generate_tokens(force=force)

    def create_db(self, force=False):
        tsv_path = os.path.join(self.folder, f"{self.name}.list.tsv")
        with open(tsv_path, "w") as f:
            iterable = self.list_dir_or_dataset
            if type(self.list_dir_or_dataset) is list:
                iterable = self.list_dir_or_dataset
            elif type(self.list_dir_or_dataset) is str:
                assert os.path.exists(self.list_dir_or_dataset), f"Folder {self.list_dir_or_dataset} does not exist"
                assert os.path.isdir(self.list_dir_or_dataset), f"{self.list_dir_or_dataset} is not a folder"
                iterable = os.listdir(self.list_dir_or_dataset)
            elif isinstance(self.list_dir_or_dataset, StructureDataset):
                iterable = self.list_dir_or_dataset.paths()

            for path in iterable:
                line = path+"\n"
                f.write(line)

        cmd = [self.foldseek_command, "createdb", tsv_path, self.db_path, "-v", self.verbose]
        print(" ".join(cmd))
        subprocess.run(cmd, check=True)
        return self

    def generate_tokens(self, force=False):
        if force or self.token_fasta_path is None:
            cmd = [self.foldseek_command, "lndb", self.db_path+"_h", self.db_path+"_ss_h", "-v", self.verbose]
            print(" ".join(cmd))
            subprocess.run(cmd , check=True)
            #cmd = [self.foldseek_command, "createindex", self.db_path, os.path.join(TEMP_FOLDER, "foldseekk"), "-v", self.verbose]
            #print(" ".join(cmd))
            #subprocess.run(cmd)
            fasta_path = self.db_path+".tokens.fasta"
            cmd = [self.foldseek_command, "convert2fasta", self.db_path+"_ss", fasta_path, "-v", self.verbose]
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)
            self.token_fasta_path = fasta_path
        return self

    def tokens_fasta(self, force=False):
        if self.token_fasta_path is None or force:
            self.generate_tokens(force=force)
        return FASTA(self.token_fasta_path)

    def sequences_fasta(self, force=False):
        if force or self.sequence_fasta_path is None:
            fasta_path = self.db_path+".aa.fasta"

            cmd = [self.foldseek_command, "convert2fasta", self.db_path, fasta_path, "-v", self.verbose]
            print(" ".join(cmd))
            subprocess.run(cmd, check=True)
            self.sequence_fasta_path = fasta_path
        
        return FASTA(self.sequence_fasta_path)

    def __iter__(self):
        self._cached_token_fasta = self.tokens_fasta()
        self._cached_aa_fasta = self.sequences_fasta()
        assert len(self._cached_token_fasta) == len(self._cached_aa_fasta), f"Token entries ({len(self._cached_token_fasta)}) and sequence entries ({len(self._cached_aa_fasta)}) do not match"
        self.i = 0
        self.max_i = len(self._cached_token_fasta)
        return self

    def __next__(self):
        if self.i < self.max_i:
            tok_name, tok_seq = self._cached_token_fasta.at(self.i)
            tok_seq = tok_seq[0]
            aa_name, aa_seq = self._cached_aa_fasta.at(self.i)
            aa_seq = aa_seq[0]
        else:
            raise StopIteration()
        self.i += 1
        assert aa_name == tok_name, f"Sequence id ({aa_name}) and Token id ({tok_name}) do not match"
        assert len(aa_seq) == len(tok_seq), f"Sequence len ({len(aa_seq)}) and Token len ({len(tok_seq)}) do not match"

        return {
            "name": aa_name,
            "tok_seq": tok_seq,
            "aa_seq": aa_seq,
        }

    @staticmethod
    def parse_name(full_name, as_dict=False):
        name = full_name.split(" ")[0]
        description = " ".join([n for n in full_name.split(" ")[1:]])
        components = name.split("_")
        chain = "-"
        model = "1"
        if len(components) == 1:
            code = components[0]
        elif len(components) == 2:
            code, chain = components
        elif (len(components) == 4) and components[1] == "MODEL":
            code, _, model, chain = components
        else:
            raise NameParsingError(f"Unable to fetch name and code from {full_name}")
        if as_dict:
            return dict(
                full_name=full_name,
                name=name,
                description=description,
                code=code,
                chain=chain,
                model=model,
            )
        else:
            return code, chain, model

    def match_dataset(self, dataset, saprot=True, atoms=True, entity_class=BIEntity, **kwargs):


        iterables = [self.tokens_fasta(**kwargs), self.sequences_fasta(**kwargs)]
        if saprot:
            iterables.append(self.saprot_embeddings(**kwargs))

        entity = None
        for (name, t_seq), (_, aa_seq), *data in zip(*iterables):

            #print(name)
            code, chain, model = self.parse_name(name)

            #print(code, chain)
            t_seq = t_seq[0]
            aa_seq = aa_seq[0]
            print(code, chain, len(t_seq), len(aa_seq))

            entry = dataset.get(code)
            try:

                #print(entry)

                try:
                    assert entity is not None
                    assert entity.code() == code
                except:
                    entity = entity_class.from_file(entry["path"], code=code, no_atoms=not atoms)


                entity.set_model(model)
                print(entity)
                print(chain)
                chains = entity.chains(chain, by_complex=True)
                #print(chains)

                sequence = ""
                for chain_entity in chains:
                    print(chain, chain_entity, chain_entity.id(), chain_entity.complex(), len(chain_entity.residues()))
                    sequence += "".join([ch.sequence() for ch in chains])
                    pass


                print(len([ch.sequence() for ch in chains]), len(sequence), len(t_seq))
                print()


            except (StructureLoadException, AssertionError) as e:
                dataset.add_to_blacklist(entry["path"], error=e)
                continue


            yield {
                "code": code,
                "chain": chain,
                "model": model,
                "t_seq": t_seq,
                "error": False,
                "aa_seq": sequence,
                "match_len": len(sequence) == len(t_seq),
                "chain_id": chain,
                "chains": chains,
                "entity": entity,
            }, *data

    @staticmethod
    def zip(aas, toks, tokens_first=False):
        assert len(aas) == len(toks), f"AAs ({len(aas)}) and Tokens ({len(toks)}) length do not match"
        l = []
        z = zip(aas, toks)
        for a, t in z:
            if tokens_first:
                l.append(f"{t.lower()}{a.upper()}")
            else:
                l.append(f"{a.upper()}{t.lower()}")
        return l


    def saprot_embeddings(self, sequence_only=False, saprot_model=None, force=False, save_folder=None, return_tensor=True):

        log(2, "Generating SaProt Embeddings...")
        if saprot_model is None:
            saprot_model = self.saprot_model
            saprot_path = f"westlake-repl/{self.saprot_model}"
        else:
            saprot_path = saprot_model
            saprot_model = saprot_path.split("/")[-1]
        from transformers import EsmTokenizer, EsmForMaskedLM
        import torch
        from src.bioiain.machine import DEVICE




        tokenizer_path = os.path.join(SUBDIR_NAME, "hf", "SaProt", f"tok_{saprot_model}")
        if not os.path.exists(tokenizer_path):
            log(3, "Downloading tokeniser:", saprot_path)
            tokenizer = EsmTokenizer.from_pretrained(saprot_path)
            os.makedirs(os.path.dirname(tokenizer_path), exist_ok=True)
            tokenizer.save_pretrained(tokenizer_path)
        tokenizer = EsmTokenizer.from_pretrained(tokenizer_path)

        model_path = os.path.join(SUBDIR_NAME, "hf", "SaProt", f"mod_{saprot_model}")
        if not os.path.exists(model_path):
            log(3, "Downloading model:", saprot_path)
            model = EsmForMaskedLM.from_pretrained(saprot_path)
            os.makedirs(os.path.dirname(model_path), exist_ok=True)
            model.save_pretrained(model_path)
        model = EsmForMaskedLM.from_pretrained(model_path)
        
        model.eval()
        model.to(DEVICE)


        for entry in self:
            log("header", entry["name"])
            if save_folder is None:
                save_folder = os.path.join(SUBDIR_NAME, "embeddings")
            #print(save_folder, "saprot", saprot_model)
            save_path = os.path.join(save_folder, "SaProtRaw", saprot_model)
            os.makedirs(save_path, exist_ok=True)
            save_path = os.path.join(save_path, entry["name"].split(" ")[0]+".pt")
            if (not force) and os.path.exists(save_path):
                if return_tensor:
                    yield torch.load(save_path), save_path, entry,  saprot_model,
                else: 
                    yield save_path, entry,  saprot_model,

                continue
            log(2, "Running SaProt model...")
            if sequence_only:
                seq = "".join([f"{aa.upper()}#" for aa in entry["aa_seq"]])
            else:
                seq = "".join(self.zip(entry["aa_seq"], entry["tok_seq"]))
            log(3, "LEN SEQ:", len(entry["aa_seq"]))
            #tokens = tokenizer.tokenize(seq)
            inputs = tokenizer(seq, return_tensors="pt")
            inputs = {k: v.to(DEVICE) for k, v in inputs.items()}
            with torch.no_grad():
                log(3, "Generating output...")
                outputs = model(**inputs, output_hidden_states=True)
                log(3, "Output ready")
                last_hidden = outputs.hidden_states[-1][:,1:-1,:]
                log(3, "EMBEDDING:", last_hidden.shape)

                torch.save(last_hidden, save_path)
                if return_tensor:
                    yield last_hidden, save_path, entry, saprot_model
                else:
                    yield save_path, entry, saprot_model

                continue


    def prost5_derive(self, from_tokens=False):
        from transformers import T5Tokenizer, T5EncoderModel
        from src.bioiain.machine import DEVICE
        import torch
        import re

        tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)

        #print("TOKENISER", tokenizer)




        seqs = []
        for entry in self:
            print(entry)
            if from_tokens:
                seqs.append(f"<fold2AA> {''.join([s.lower() for s in entry['tok_seq']])}")
            else:
                seqs.append(f"<AA2fold> {''.join([s.upper() for s in entry['aa_seq']])}")

            ids = tokenizer.batch_encode_plus(seqs,
                add_special_tokens=True,
                padding="longest",
                return_tensors='pt').to(DEVICE)
            print(ids.input_ids)


        model = T5EncoderModel.from_pretrained("Rostlab/ProstT5").to(DEVICE)
        print("MODEL", model)

        model.float() if DEVICE=='cpu' else model.half()
        with torch.no_grad():
            embedding_repr = model(
                ids.input_ids, 
                attention_mask=ids.attention_mask
                )
            for n in range(embedding_repr.last_hidden_state.shape[0]):
                emb = embedding_repr.last_hidden_state[n]
                print("Embedding", n)
                print(emb)















