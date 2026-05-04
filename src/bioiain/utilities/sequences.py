import os, sys, json, subprocess
import tempfile

from ..utilities.exceptions import *
from .logging import log
from .. import TEMP_FOLDER, SUBDIR_NAME


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
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.single_line = None


    def __repr__(self):
        return f"<bi.{self.__class__.__name__}: {self.fasta_path}>"


    def _parse_fasta(self, names=True, sequences=True, key=None):
        assert names or sequences

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
    def __init__(self, *args, mmseqs_cmd="mmseqs", db_name=None, verbosity=2, **kwargs):
        super().__init__(*args, **kwargs)
        self.tmp_folder = os.path.join(TEMP_FOLDER, "mmseqs2")
        os.makedirs(self.tmp_folder, exist_ok=True)
        self.mmseqs_cmd = mmseqs_cmd
        self.db_name = None
        self.db_folder = None
        self.verbosity = verbosity
        self.name = self.name.replace(".fasta", ".db")
        if db_name is None:
            db_name = self.name.split(".")[0]
        if os.path.isdir(self.fasta_path):
            log(2, "Input is already DB, setup only")
            self.db_folder = self.fasta_path
            self.db_name = db_name
            self.db_path = os.path.join(self.db_folder, self.db_name)
        else:
            log(2, "Input is a file, creating DB...")
            self.create_db(db_name=db_name, **kwargs)



    def _cmd(self, command, *args, **kwargs):

        cmd = [self.mmseqs_cmd, command]

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




    def create_db(self, db_name=None, fasta_path=None, **kwargs):
        if fasta_path is None:
            self.fasta.rewrite(key_start=">")
            fasta_path = self.fasta_path
        if db_name is None:
            db_name = self.name.split(".")[0]
        db_folder = os.path.join(SUBDIR_NAME, "mmseqs", self.name)
        db_path = os.path.join(db_folder, db_name)
        os.makedirs(db_folder, exist_ok=True)
        self._cmd("createdb", fasta_path, db_path, createdb_mode=0)
        self.db_name = db_name
        self.db_folder = db_folder
        self.db_path = db_path
        return self


    def cluster(self, db_name=None, reassign=False, force=False, linear=False, easy=False, **kwargs):
        if db_name is None:
            db_name = self.db_name
        cluster_db_folder = os.path.join(self.db_folder.replace(".db", ".cluster"))
        cluster_db_path = os.path.join(cluster_db_folder, db_name)
        out_path = os.path.join(cluster_db_folder, f"{db_name}_clustered.tsv")
        data_path = out_path.replace(".tsv", ".json")
        if linear:
            cmd = ["linclust"]
        else:
            cmd = ["cluster"]
        if easy:
            cmd =  ["easy-"+cmd[0]]
        cmd.extend([self.db_path, cluster_db_path, self.tmp_folder])
        if reassign and not linear:
            cmd.append("--cluster-reassign")

        params = {
            "cmd": " ".join([str(c) for c in cmd]),
            "reassign":reassign,
            "linear":linear,
            "easy":easy,
        }
        os.makedirs(cluster_db_folder, exist_ok=True)

        if not os.path.exists(cluster_db_path) or not os.path.exists(data_path):
            force=True
        if os.path.exists(data_path):
            if json.load(open(data_path))["params"] != params:
                log(3, "Different params detected")
                force = True

        if force:
            try:
                self._cmd(*cmd, v=self.verbosity)
            except:
                raise ClusteringError()
        else:
            log(3, "Cluster DB already clustered (mmseqs2)")

        if force or not os.path.exists(out_path):
            try:
                self._cmd("createtsv", self.db_path, self.db_path, cluster_db_path, out_path, v=self.verbosity)
            except:
                raise ClusteringError()

        try:
            clusters = {}
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
        return data_path








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









