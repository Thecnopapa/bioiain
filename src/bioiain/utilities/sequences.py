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



    def _parse_fasta(self, names=True, sequences=True):
        assert names or sequences

        fasta_dict = {}
        with open(self.fasta_path) as f:
            next_seq = False
            last_key = None
            for line in f.readlines():
                line = line.replace("\n", "").strip()
                if line.startswith("#"):
                    next_seq = True
                    continue
                if line.startswith(">"):
                    name = line[1:]
                    if name not in fasta_dict:
                        fasta_dict[name] = []
                    next_seq = True
                    last_key = name
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
                        fasta_dict[last_key] += line


        #print(names, sequences)
        if names and sequences:
            return fasta_dict
        elif names:
            return list(fasta_dict.keys())
        elif sequences:
            seqs = []
            [seqs.extend(seq) for seq in fasta_dict.values()]
            return seqs





    def get_names(self):
        return self._parse_fasta(names=True, sequences=False)
    def get_sequences(self):
        return self._parse_fasta(names=False, sequences=True)
    def parse(self):
        return self._parse_fasta()





class MSA(object):
    def __init__(self, fasta_path, name=None):
        self.fasta_path = fasta_path
        fasta = FASTA(fasta_path)
        self.fasta_dict = fasta.parse()
        if name is None:
            name = os.path.basename(fasta_path)
        self.name = name
        log("header", f"Initialising MSA: {self}")
        self.msa_path = self._run_clustal_msa(name=self.name)
        self.tree_path = self._build_tree(self.msa_path)

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}:{self.name} ({len(self)} sequences)>"

    def __len__(self):
        return len(self.fasta_dict)



    def _run_clustal_msa(self, fasta_path=None, name="temp", out_folder=None, clustal_cmd="clustalw", matrix="BLOSUM", out_format="fasta"):
        
        if fasta_path is None:
            fasta_path = self.fasta_path
        log(2, f"Calculating MSA of: {fasta_path}")
        fname = f"{name}_{matrix}.ms.alignment.fasta"
        if out_folder is None:
            out_folder = "/tmp/bioiain/alignments"
        os.makedirs(out_folder, exist_ok=True)
        out_path = os.path.join(out_folder, fname)
        cmd = [
            "clustalw", "-align", "-type=protein",
            f"-infile={fasta_path}",
            f"-matrix={matrix}",
            f"-outfile={out_path}",
            f"-output={out_format}"
        ]

        #print("$", " ".join(cmd))
        out_log = open("/dev/null", "w")
        subprocess.run(cmd, stdout=out_log)
        return out_path
        

    def _build_tree(self, align_path):
        log(2, f"Building tree for: {align_path}")
        cmd = [
            "clustalw", "-tree", "-type=protein",
            f"-infile={align_path}",
            "-outputtree=nj",
        ]

        #print("$", " ".join(cmd))
        out_path = align_path.replace(".fasta", ".nj")
        comp_file = out_path + ".list"

        f = open(comp_file, "w")
        subprocess.run(cmd, stdout=f)
        return out_path




    def get_similar(self, target, name="temp", similarity=95):
        threshold = (100-similarity) / 100
        log(2, f"Finding similar at {similarity}% for {target}")

        seq_num = self._get_seq_num(target)
        neighbour_nums = self._get_neighbours(seq_num)
        neighbour_names = [self._get_seq_name(n) for n in neighbour_nums]

    
        exit()
        print(neighbour_names)
        return neighbour_names


    def _get_seq_num(seq_name):
        log(2, f"Finding seq_num for {seq_name}")
        comp_path = self.tree_path+".list"
        seq_num = None

        #TODO: Parse seq nums
        return seq_num

    def _get_seq_name(seq_num):
        log(2, f"Finding seq_name for {seq_num}")
        comp_path = self.tree_path+".list"
        seq_name = None

        #TODO: Parse seq names
        return seq_name


    def _get_neighbours(seq_num, threshold=0.05):
        log(2, f"Finding neighbours (seq. {seq_num}), threshold={threshold}")
        neighbours = []

        #TODO: Parse neighbours
        return neighbours
        








             