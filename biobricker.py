seq = "ATGAAACCC"
seq2 = "GAATTCGCGGCCGCTTCTAGATGAAACCCTACTAGTAGCGGCCGCTGCAG"
seq3 = "ALASDAIR"
seq4 = "GAATTCGCGGCCGCTTCTAGAGALASDAIRTACTAGTAGCGGCCGCTGCAG"

class Sequence:
    def __init__(self, name, sequence):
        self.name = name
        self.sequence = sequence
    def biobrick(self):
        if self.sequence[:3]=="ATG": #identifies coding sequences
            bb = "GAATTCGCGGCCGCTTCTAG" + self.sequence + "TACTAGTAGCGGCCGCTGCAG"
        else:
            bb = "GAATTCGCGGCCGCTTCTAGAG" + self.sequence + "TACTAGTAGCGGCCGCTGCAG"
        return bb
    def debiobrick(self):
        if self.sequence[:20] == "GAATTCGCGGCCGCTTCTAG":
            if self.sequence[:22] == "GAATTCGCGGCCGCTTCTAGAG":
                dbb = self.sequence[22:-21]
            else:
                dbb = self.sequence[20:-21]
            return dbb
        else: print "Not a BioBrick"
        
seq1=Sequence("seq1",seq)
seq2=Sequence("seq2",seq2)
seq3=Sequence("seq2",seq3)
seq4=Sequence("seq2",seq4)