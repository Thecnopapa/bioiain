import os, sys, json, subprocess


from .logging import log


d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}








class FASTA(object):
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path


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





    def get_names(self, key=None):
        return self._parse_fasta(names=True, sequences=False, key=key)

    def get_sequences(self, key=None):
        return self._parse_fasta(names=False, sequences=True, key=key)

    def parse(self, key=None):
        return self._parse_fasta(key=key)






class MSA(object):
    def __init__(self, fasta_path, name=None, verbose=False):
        self.fasta_path = fasta_path
        fasta = FASTA(fasta_path)
        self.fasta_dict = fasta.parse()
        if name is None:
            name = os.path.basename(fasta_path)
        self.name = name
        log("header", f"Initialising MSA: {self}")
        self.msa_path = self._run_clustal_msa(name=self.name, verbose=verbose)
        self.tree_path = self._build_tree(self.msa_path)
        self.msa_fasta = FASTA(self.msa_path)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} ({len(self)} sequences)>"

    def __len__(self):
        return len(self.fasta_dict)



    def _run_clustal_msa(self, fasta_path=None, name="temp", out_folder=None, clustal_cmd="clustalw", matrix="BLOSUM", out_format="fasta", force=False, verbose=False):

        if fasta_path is None:
            fasta_path = self.fasta_path
        log(2, f"Calculating MSA of: {fasta_path}")
        fname = f"{name}_{matrix}.ms.alignment.fasta"
        if out_folder is None:
            out_folder = "/tmp/bioiain/alignments"
        os.makedirs(out_folder, exist_ok=True)
        out_path = os.path.join(out_folder, fname)
        if os.path.exists(out_path) and not force:
            log(3, "Alignment already generated")
            return out_path
        cmd = [
            "clustalw", "-align", "-type=protein",
            f"-infile={fasta_path}",
            f"-matrix={matrix}",
            f"-outfile={out_path}",
            f"-output={out_format}"
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
                        comps = [l for l in re.split(" |.|;|=",line.strip()) if l != ""]
                        num1 = int(comps[0])
                        num2 = int(comps[2])
                        dist = float(comps[5].replace(";", ""))
                        length = int(comps[8].replace("\n", ""))
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









