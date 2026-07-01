import os, sys, json, subprocess, shutil

from ..utilities.exceptions import *
from .logging import log
from .. import TEMP_FOLDER, SUBDIR_NAME
import polars as pl


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', "UNK": "X"}

d3toint =  {'CYS': 0, 'ASP': 1, 'SER': 2, 'GLN': 3, 'LYS': 4,
             'ILE': 5, 'PRO': 6, 'THR': 7, 'PHE': 8, 'ASN': 9,
             'GLY': 10, 'HIS': 11, 'LEU': 12, 'ARG': 13, 'TRP': 14,
             'ALA': 15, 'VAL': 16, 'GLU': 17, 'TYR': 18, 'MET': 19, "UNK":20}

def d3(resname):
    try:
        ri = d3toint[resname]
        rn = d3to1[resname]
    except:
        log("warning", f"Unknown resname: {resname} (using UNK/X)")
        ri = 20
        rn = "X"
    return rn, ri



def intto1(i):
    for d3, n in d3toint.items():
        if i == n:
            return d3to1[d3]
    return "X"




class FASTA(object):
    def __init__(self, fasta_path, force=False, use_cache=True):
        self.fasta_path = fasta_path
        self.single_line = None
        self._cache = None
        self.force = force
        self.use_cahce = use_cache


    def __repr__(self):
        return f"<bi.{self.__class__.__name__}: {self.fasta_path}>"


    def _parse_fasta(self, names=True, sequences=True, key=None, use_cache=None, force=None):
        assert names or sequences

        if force is None:
            force = self.force
            use_cache = False
        if use_cache is None:
            use_cache = self.use_cache

        if use_cache and self._cache is not None and key is None:
            fasta_dict = self._cache
        else:
            if key is not None:
                if type(key) is str:
                    key = [key]
                elif type is not list:
                    key = list(key)

            fasta_dict = {}
            with open(self.fasta_path) as f:
                next_seq = False
                last_key = None
                wait_key = False
                for line in f.readlines():
                    line = line.replace("\n", "").strip()
                    if line.startswith("#"):
                        next_seq = True
                        continue
                    if line.startswith(">"):
                        wait_key = False
                        name = line[1:].strip()
                        if key is not None:
                            if len(key) == 0:
                                break
                            #print(name, key, name in key)
                            if name in key:
                                key.remove(name)
                            else:
                                wait_key = True
                                continue
                        if name not in fasta_dict:
                            fasta_dict[name] = []
                        next_seq = True
                        last_key = name
                        continue
                    elif wait_key:
                        continue
                    if not sequences:
                        continue
                    if line.strip() == "":
                        next_seq = True
                        continue
                    else:
                        if last_key is None:
                            continue
                        if next_seq:
                            fasta_dict[last_key].append(line)
                            next_seq = False
                        else:
                            fasta_dict[last_key][-1] += line

            if key is not None:
                #print(key)
                assert len(key) == 0

            if use_cache:
                self._cahche = fasta_dict

        #print(names, sequences)
        if names and sequences:
            return fasta_dict
        elif names:
            return list(fasta_dict.keys())
        elif sequences:
            seqs = []
            [seqs.extend(seq) for seq in fasta_dict.values()]
            return seqs

    def rewrite(self, duplicates=False, empties=False, space_between=False, key_start="> "):
        log(3, "Rewriting FASTA:", self.fasta_path)
        data = self._parse_fasta()
        with open(self.fasta_path, "w") as f:
            for key, sequences in data.items():
                if not empties:
                    sequences = [s for s in sequences if len(s) > 0]
                n_seqs = len(sequences)
                if not duplicates:
                    if n_seqs > 1:
                        log("warning", f"{n_seqs} sequences for id: {key} (keeping only first)")
                    sequences = sequences[:1]

                for seq in sequences:
                    f.write(f"{key_start}{key}\n")
                    f.write(f"{seq}\n")
                    if space_between:
                        f.write("\n")
        self.single_line = True
        return self.fasta_path

    def get_names(self, key=None):
        return self._parse_fasta(names=True, sequences=False, key=key)

    def get_sequences(self, key=None):
        return self._parse_fasta(names=False, sequences=True, key=key)


    def parse(self, key=None):
        return self._parse_fasta(key=key)

    def keys(self):
        return self.get_names()

    def values(self):
        return self.get_sequences()

    def sequences():
        return self.get_sequences()

    def dict():
        return self.parse()

    def items():
        return self.dict().items()

    def index(self, key):
        return self.keys().index(key)

    def at(self, index):
        key = self.keys()[index]
        return key, self.get_sequences(key)

    def __len__(self):
        return len(self.keys())

    def __iter__(self):
        self.i = 0
        return self

    def __next__(self):
        if self.i < len(self.keys()):
            self.i += 1
            return self.at(self.i-1)
        else:
            raise StopIteration()





class MSA(object):
    def __init__(self, fasta_path, name=None, **kwargs):
        self.fasta_path = fasta_path
        self.fasta = FASTA(fasta_path)
        if name is None:
            name = os.path.basename(fasta_path).replace(".fasta", "")
        self.name = name
        log(1, f"Initialising {self.__class__.__name__}...")
        log(2, "Fasta path:", self.fasta_path)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} ({len(self)} sequences)>"


    def __len__(self):
        return len(self.fasta.get_names())





class MMSEQS2(MSA):
    def __init__(self, *args, mmseqs_cmd="mmseqs", db_name=None, verbosity=1, folder=None, force=False, **kwargs):
        super().__init__(*args, **kwargs)
        self.fasta.rewrite(key_start=">")
        self.tmp_folder = os.path.join(TEMP_FOLDER, "mmseqs2")
        os.makedirs(self.tmp_folder, exist_ok=True)
        self.databases = {}
        self.mmseqs_cmd = mmseqs_cmd

        self.verbosity = verbosity
        self.name = self.name.replace(".dataset", "")
        if not self.name.endswith(".mmseqs"):
            self.name += ".mmseqs"

        if folder is None:
            folder = os.path.join(SUBDIR_NAME, "mmseqs")
        self.db_folder = os.path.join(folder, self.name)

        os.makedirs(self.db_folder, exist_ok=True)

        if force:
            self.delete()

        if db_name is None:
            db_name = self.name.split(".")[0]
        self.db_name = db_name

        if os.path.exists(self.db_path()):
            log(2, "Input is already DB, setup only")
        else:
            log(2, "Input is a file, creating DB...")
            self.create_db(db_name=db_name, **kwargs)

    def db_path(self, suffix="db"):
        return os.path.join(self.db_folder, f"{self.db_name}.{suffix}")

    def delete(self, db_name=None, make_dir=True):
        if db_name is None:
            shutil.rmtree(self.db_folder, ignore_errors=True)
        else:
            for file in self.db_folder:
                prefix = self.db_path(suffix=db_name)
                if file.startswith(prefix):
                    os.remove(os.path.join(self.db_folder, file))
        os.makedirs(self.db_folder, exist_ok=True)


    def _cmd(self, command, *args, **kwargs):

        if type(command) is str:
            command = [command]

        cmd = [self.mmseqs_cmd, *command]
        if "v" not in kwargs:
            kwargs["v"] = kwargs.pop("verbosity", self.verbosity)

        for kwarg, value in kwargs.items():
            if not kwarg.startswith("--"):
                if len(kwarg) == 1:
                    kwarg = f"-{kwarg}"
                else:
                    kwarg = f"--{kwarg.replace('_', '-')}"
            cmd.extend([kwarg, str(value)])

        cmd.extend([str(a) for a in args])
        log(3, "$", " ".join(cmd))
        subprocess.run(cmd)




    def create_db(self, force=False, fasta_path=None, **kwargs):
        if force:
            self.delete()
        if fasta_path is None:
            self.fasta.rewrite(key_start=">")
            fasta_path = self.fasta_path

        self.databases["sequence"] = self.db_path()
        self._cmd("createdb", fasta_path, self.databases["sequence"], createdb_mode=1, shuffle=0)
        return self


    def  write(self, query, target=None, result=None, output_file=None, mode="tab", **kwargs):

        if mode.lower() == "tsv":
            cmd = ["createtsv"]
            extension = "tsv"
        elif mode.lower() == "fasta":
            cmd = ["result2flat"]
            extension = "fasta"
        elif mode.lower() == "tab" or mode.lower() == "alis":
            cmd = ["convertalis"]
            extension = "tab"
        else:
            raise NotImplementedError


        cmd.append(query)
        if output_file is None:
            if not output_file.endswith(f".{extension}"):
                output_file += f".{extension}"
            output_file = ".".join(query.split(".")[:-1])
        if target is not None:
            cmd.append(target)
        if result is not None:
            cmd.append(result)
        cmd.append(output_file)

        try:
            self._cmd(cmd, **kwargs)
        except:
            raise TsvError()
        return output_file

    def cluster(self, reassign=False, force=False, linear=False, easy=False, **kwargs):


        self.databases["clustered"] = self.db_path("cluster")
        out_path = self.db_path("cluster")+ ".tsv"
        fasta_path = self.db_path("cluster")+ ".fasta"
        data_path = self.db_path("cluster") + ".json"

        if linear:
            cmd = ["linclust"]
        else:
            cmd = ["cluster"]
        if easy:
            cmd =  ["easy-"+cmd[0]]

        cmd.extend([self.db_path(), self.db_path("cluster"), self.tmp_folder])

        if reassign and not linear:
            cmd.append("--cluster-reassign")

        params = {
            "cmd": " ".join([str(c) for c in cmd]),
            "reassign":reassign,
            "linear":linear,
            "easy":easy,
        }

        if not os.path.exists(self.db_path("cluster")) or not os.path.exists(data_path):
            force=True
        if os.path.exists(data_path):
            if json.load(open(data_path))["params"] != params:
                log(3, "Different params detected")
                force = True

        if force:
            self.delete("cluster")
            try:
                self._cmd(*cmd)
            except:
                raise ClusteringError()
        else:
            log(3, "Cluster DB already clustered (mmseqs2)")

        if force or not os.path.exists(out_path):
            try:
                tsv = self.write(self.db_path(), self.db_path(), self.db_path("cluster"), out_path, mode="tsv")

            except:
                raise TsvError()

        try:
            clusters = {}
            print(out_path)
            with open(out_path) as f:
                for line in f:
                    c, i = line.strip().split("\t")
                    if c not in clusters:
                        clusters[c] = {"name":c, "list": []}
                    clusters[c]["list"].append(i)
            data = {"params": params, "clusters":{},}
            for c in clusters:
                data["clusters"][len(data["clusters"])] = {**clusters[c], "n": len(clusters[c]["list"])}

            json.dump(data, open(data_path, "w"), indent=4)
        except:
            raise ClusteringError()
        print("Cluster data:", data_path)
        return data_path

    def map(self, *args, **kwargs) -> pl.DataFrame:
        kwargs.pop("map", None)
        return self.search(*args, map=True, **kwargs)

    def search(self, query_db, exhaustive=True, map=False, dataset_name=None, **kwargs) -> pl.DataFrame:

        if map:
            cmd = ["map"]
            folder_name = "map"
        else:
            cmd = ["search"]
            folder_name = "search"
            if exhaustive:
                cmd.append("--exhaustive-search")
                cmd.extend(["--alignment-mode", "3"])

        log(1, f"Searching({query_db}) in {self.db_path()} cmd={folder_name}")



        aligned_db = os.path.join(self.tmp_folder, folder_name, str(dataset_name))
        shutil.rmtree(aligned_db, ignore_errors=True)
        os.makedirs(aligned_db, exist_ok=True)
        aligned_db = os.path.join(aligned_db, f"temp.{folder_name}")



        cmd.extend([query_db, self.db_path(), aligned_db, self.tmp_folder])




        cmd.append("-a")

        try:
            self._cmd(*cmd, **kwargs)
        except:
            raise SearchError()

        columns = "query,target,evalue,raw,bits,fident,alnlen,pident,qcov,tcov,qlen,tlen,qstart,tstart,qaln,taln"
        tsv = self.write(query_db, self.db_path(), aligned_db, output_file=".".join(query_db.split(".")[:-1]) + f".{folder_name}.tab", mode="alis", format_mode=4, format_output=columns)
        df = pl.read_csv(tsv, separator="\t", schema_overrides={"query":pl.String, "target":pl.String})
        log(df)
        return df





class CLUSTAL(MSA):
    def __init__(self, *args, verbose=False, run_msa=True, build_tree=False, matrix_path=None, **kwargs):
        super().__init__(*args, **kwargs)
        kwargs.pop("name", None)
        if run_msa:
            self.msa_path = self._run_clustal_msa(name=self.name, verbose=verbose, matrix_path=matrix_path, **kwargs)
            self.msa_fasta = FASTA(self.msa_path)
            self.msa_fasta.rewrite()
        if build_tree:
            self.tree_path = self._build_tree(self.msa_path)


    def _run_clustal_msa(self, fasta_path=None, name="temp", out_folder=None, clustal_cmd="clustalw", matrix="BLOSUM", out_format="fasta", force=False, verbose=False, matrix_path=None, **kwargs):

        if fasta_path is None:
            fasta_path = self.fasta_path
        log(2, f"Calculating MSA ({matrix}) of: {fasta_path}")
        fname = f"{name}_{matrix}.ms.alignment.fasta"
        if matrix == "path":
            assert matrix_path is not None
            matrix = matrix_path
        if out_folder is None:
            out_folder = os.path.join(TEMP_FOLDER, "alignments")
        os.makedirs(out_folder, exist_ok=True)
        out_path = os.path.join(out_folder, fname)
        if os.path.exists(out_path) and not force:
            log(3, "Alignment already generated (CLUSTAL)")
            return out_path
        cmd = [
            clustal_cmd, "-align", "-type=protein",
            f"-infile={fasta_path}",
            f"-matrix={matrix}",
            f"-pwmatrix={matrix}",
            f"-outfile={out_path}",
            f"-output={out_format}"
            f"-slow"
        ]
        if verbose:
            log(3, "$", " ".join(cmd))
            subprocess.run(cmd)
        else:
            out_log = open("/dev/null", "w")
            subprocess.run(cmd, stdout=out_log)
        return out_path


    def _build_tree(self, align_path, force=False):
        log(2, f"Building tree for: {align_path}")
        out_path = align_path.replace(".fasta", ".nj")
        comp_file = out_path + ".list"
        if os.path.exists(out_path) and os.path.exists(comp_file) and not force:
            log(3, "Tree already generated")

            return out_path
        cmd = [
            "clustalw", "-tree", "-type=protein",
            f"-infile={align_path}",
            "-outputtree=nj",
        ]

        #print("$", " ".join(cmd))


        f = open(comp_file, "w")
        subprocess.run(cmd, stdout=f)
        return out_path


    def get_similar(self, target, name="temp", similarity=95):
        threshold = (100-similarity) / 100
        log(2, f"Finding similar at {similarity}% for {target}")

        seq_num = self._get_seq_num(target)
        neighbour_nums = self._get_neighbours(seq_num, threshold=threshold)
        neighbour_names = [self._get_seq_name(n) for n in neighbour_nums]

        #print(neighbour_names)

        #exit()
        return neighbour_names


    def _get_seq_num(self, seq_name) -> int|None:
        #log(3, f"Finding seq_num for {seq_name}")
        comp_path = self.tree_path+".list"
        seq_num = None
        with open(comp_path, "r") as f:
            for line in f.readlines():
                comps = line.split(" ")
                if len(comps) < 2:
                    continue
                if seq_name in comps:
                    seq_num = int(comps[1].replace(":", ""))
                    break
        return seq_num


    def _get_seq_name(self, seq_num):
        #log(3, f"Finding seq_name for {seq_num}")
        comp_path = self.tree_path+".list"
        seq_name = None
        with open(comp_path, "r") as f:
            for line in f.readlines():
                comps = line.split(" ")
                if len(comps) < 2:
                    continue
                if f"{seq_num}:" == comps[1]:
                    seq_name = comps[2]
        return seq_name


    def _get_neighbours(self, seq_num, threshold=0.05):
        log(3, f"Finding neighbours (seq. {seq_num}), threshold={threshold}")
        import re
        neighbours = []
        with open(self.tree_path, "r") as f:
            for line in f.readlines():
                if "DIST" in line and "length" in line:
                    try:
                        comps = [l for l in re.split(' |vs\.|;|=', line.strip()) if l != ""]
                        num1 = int(comps[0])
                        num2 = int(comps[1])
                        dist = float(comps[3])
                        length = int(comps[5].replace("\n", ""))
                        if dist > threshold:
                            continue
                        if seq_num == num1:
                            neighbours.append(num2)
                        elif seq_num == num2:
                            neighbours.append(num1)
                    except Exception as e:
                        log("warning", f"Error reading tree file: {self.tree_path}")
                        print(line)
                        print(comps)
                        raise e



        log(3, f"Found {len(neighbours)} neighbours")
        return neighbours








