import os, sys, json





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
    def __init__(self, fasta_path):
        self.fatsa_path = fasta_path
        fasta = FASTA(fasta_path)
        self.fasta_dict = fasta.parse()

    def __repr__(self):
        return f"<bi.{self.__class__.__name__}: {len(self)} sequences>"

    def __len__(self):
        return len(self.fasta_dict)


    def get_similar(self, target, similarity=95):
        print(f"Finding similar at {similarity}%")
        print("Target:", target)
        print(self.fasta_dict)


        return self.fasta_dict


             