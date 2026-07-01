import os, sys, json, subprocess

sys.path.append('..')

from src.bioiain.base import BIEntity
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
    def __init__(self,name, list_dir_or_dataset, folder=None, foldseek_command="foldseek", force=False):
        self.foldseek_command = foldseek_command
        if folder is None:
            folder = os.path.join(SUBDIR_NAME, "foldseek")
        self.folder = os.path.join(folder, name)
        self.db_path = os.path.join(self.folder, name)
        os.makedirs(self.folder, exist_ok=True)
        self.name = name
        self.list_dir_or_dataset = list_dir_or_dataset

        self.token_fasta_path = None
        self.sequence_fasta_path = None


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

        cmd = [self.foldseek_command, "createdb", tsv_path, self.db_path, "-v", "3"]
        print(" ".join(cmd))
        subprocess.run(cmd)
        return self

    def generate_tokens(self, force=False):
        if force or self.token_fasta_path is None:
            cmd = [self.foldseek_command, "lndb", self.db_path+"_h", self.db_path+"_ss_h", "-v", "3"]
            print(" ".join(cmd))
            subprocess.run(cmd)
            #cmd = [self.foldseek_command, "createindex", self.db_path, os.path.join(TEMP_FOLDER, "foldseekk"), "-v", "3"]
            #print(" ".join(cmd))
            #subprocess.run(cmd)
            fasta_path = self.db_path+".tokens.fasta"
            cmd = [self.foldseek_command, "convert2fasta", self.db_path+"_ss", fasta_path, "-v", "3"]
            print(" ".join(cmd))
            subprocess.run(cmd)
            self.token_fasta_path = fasta_path
        return self

    def tokens_fasta(self, force=False):
        if self.token_fasta_path is None or force:
            self.generate_tokens(force=force)
        return FASTA(self.token_fasta_path)

    def sequences_fasta(self, force=False):
        if force or self.sequence_fasta_path is None:
            fasta_path = self.db_path+".aa.fasta"

            cmd = [self.foldseek_command, "convert2fasta", self.db_path, fasta_path, "-v", "3"]
            print(" ".join(cmd))
            subprocess.run(cmd)
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



    def match_dataset(self, dataset):

        fasta = self.tokens_fasta()

        entity = None
        for name, t_seq in fasta.parse().items():

            name = name.split(" ")[0]
            print(name.split("_"))
            if len(name.split("_")) == 1:
                code = name.split("_")[0]
                chain = "*"
            elif len(name.split("_")) == 2:
                code, chain = name.split("_")
            else:
                raise Exception(f"Unable to fetch name and code from {name}")
            print(code, chain)
            t_seq = t_seq[0]
            print(code, chain, len(t_seq))
                
            try:
                entry = dataset.get(code)
                print(entry)

                try:
                    assert entity is not None
                    assert entity.code() == code
                except:
                    entity = BIEntity.from_file(entry["path"], code=code)
                print(entity)

                chains = entity.chains(chain)
                for chain_entity in chains:
                    print(chain, chain_entity, chain_entity.id(), chain_entity.complex())
                assert len(chains) == 1
                ch = chains[0]

                print(len(ch.sequence()), len(t_seq))
                print()
            except:
                yield {
                    "t_seq": t_seq,
                    "error": true,
                    "chain_id": chain,
                }
            yield {
                "t_seq": t_seq,
                "error": false,
                "aa_seq": ch.sequence(),
                "match_len": len(ch.sequence()) == len(t_seq),
                "chain_id": chain,
                "chain": chain_entity,
                "entity": entity,
            }

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



    def prost5_derive(self, from_tokens=False):
        from transformers import T5Tokenizer, T5EncoderModel
        from src.bioiain.machine import DEVICE
        import torch
        import re

        tokenizer = T5Tokenizer.from_pretrained('Rostlab/ProstT5', do_lower_case=False)

        print("TOKENISER", tokenizer)




        seqs = []
        for entry in self:
            print(entry)
            if from_tokens:
                seqs.append(f"<fold2AA> {"".join([s.lower() for s in entry["tok_seq"]])}")
            else:
                seqs.append(f"<AA2fold> {"".join([s.upper() for s in entry["aa_seq"]])}")

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















class CompactStructure(FragmentedStructure):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def ca_kdtree(self, force=False, **kwargs):
        if force or getattr(self, "_kdtrees", {}).get("ca", None) is None:
            KDT(self, mode="ca", auto_parse_symmetry=True, **kwargs)
        return self._kdtrees["ca"]

    def _calculate_compactness(self, radius=10, plot=False, session=False):
        kdtree = self.ca_kdtree()
        print(kdtree)
        from src.bioiain.visualisation.plots import plasma

        if plot:
            from src.bioiain.visualisation.plots import fig3D, close, show, line
            fig, ax = fig3D()
            print(fig, ax)

        if session:
            from src.bioiain.visualisation.pymol import PymolScript
            script = PymolScript(name=f"compactness_{self.name()}", folder = self.folder())
            minimal = self.path(minimal=True)
            entity_name = script.load(minimal)
            print(script ,entity_name)


        max_frags = max([k["atom"].get_misc("fragment", 0) for k in kdtree])
        for n, k in enumerate(kdtree):
            if k["op"] != 1:
                continue
            # print(n, k)
            fragment = k["atom"].get_misc("fragment", 0)
            color = plasma(fragment, scale=max_frags)
            coord = k["coord"]
            atom = k["atom"]
            neighs = list(kdtree.radius(k["coord"], radius=radius)[0])
            # print(neighs)
            final_vector = np.array([0., 0., 0.])
            valid_nn = 0
            for nn in neighs:
                if n == nn:
                    # log("warning", "Same atom:", n, nn)
                    continue
                if (kdtree.atom_of(nn).get_misc("fragment", None) == fragment) and (kdtree.pos_of(nn) in [None, 1]):
                    # log("warning", "Same fragment:", fragment,  kdtree.atom_of(nn).get_misc("fragment", None),)
                    continue
                valid_nn += 1
                # ax.plot(*line(k["coord"], kdtree.coord_of(nn)), c=color)
                final_vector += np.array(vector(coord, kdtree.coord_of(nn)))



            if valid_nn > 0:
                final_vector /= valid_nn
                final_vector *= -1
                compactness = length(final_vector)
                # print(final_vector, compactness)
                ccol = plasma(compactness, scale=10)
                hexccol = plasma(compactness, scale=10, as_pymol_hex=True)
                vector_end = coord + final_vector
                if plot:
                    ax.scatter(*coord, c=ccol, s=valid_nn + 1)
                    ax.plot(*line(coord, vector_end), c=ccol)
                if session:
                    r1 = f"({entity_name} and i. {atom.resnum} and c. {atom.complex})"
                    #print(r1)
                    a1 = f"({r1} and n. ca)"
                    script.line(name="compactness", sele1=a1, coord2=vector_end)
                    script.color(r1, color=hexccol)


        if plot:
            show()
            close(fig)

        if session:
            script.compile()
            script.execute()













