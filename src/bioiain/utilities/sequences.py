import os, sys, json





d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
             'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
             'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
             'ALA': 'A', 'VAL': 'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}








class FASTA(object):
    def __init__(self, fasta_path):
        self.fasta_path


    def _parse_fasta(self, names=True, sequences=True):
        assert names or sequences

        fasta_dict = {}
        with open(self.fasta_path) as f:
            next_seq = False
            for line in f.read():
                line = line.replace("\n", "").strip().replace(" ","")
                if line.startswith("#"):
                    next_seq = True
                    continue
                if line.startswith(">"):
                    name = line[1:]
                    if name in fasta_dict:
                        fasta_dict[name] = []
                    next_seq = True
                    continue
                if not sequences:
                    continue
                if line.strip() == "":
                    next_seq = True
                    continue
                else:
                    if next_seq:
                        fasta_dict[fasta_dict.keys()[-1]].append(line)
                        next_seq = False
                    else:
                        fasta_dict[fasta_dict.keys()[-1]][-1] += line

        if names and sequences:
            return fasta_dict
        elif names:
            return fasta_dict.keys()
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



    def get_similar(self, target, similarity=95):
        pass


             