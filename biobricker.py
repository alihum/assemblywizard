# imports the argparse and SeqIO modules
# argparse is used to detect command line arguments
# SeqIO allows you to parse DNA files
import argparse
from Bio import SeqIO

#sets up a parser to run in the command line
parser = argparse.ArgumentParser(description="This program makes sequences into BioBricks and turns BioBricks into normal sequences")

# this allows a file name to be specified in the command line
parser.add_argument('filename',
help="Insert name of fasta file - must be in same directory")
# this argument sets an option to BioBrick a sequence
parser.add_argument("-b",action='store_true',
help="BioBrick a sequence")
# this argument sets an option to de-BioBrick a Biobrick
parser.add_argument("-d",action='store_true',
help="Produce normal sequence from a BioBrick")

# This sets up a class for a sequence object
class Sequence:

    def __init__(self,sequence):
        self.sequence = sequence
# this function takes the sequence and attaches the prefix and suffix to
# the sequence
    def biobrick(self):
# as coding sequences take a different prefix this searches for this
# and applies it if necessary
        if self.sequence[:3]=="ATG":
            bb = "GAATTCGCGGCCGCTTCTAG" + self.sequence + "TACTAGTAGCGGCCGCTGCAG"
        else:
            bb = "GAATTCGCGGCCGCTTCTAGAG" + self.sequence + "TACTAGTAGCGGCCGCTGCAG"
        return bb

# this function reverses the process of BioBricking and returns an
# unprefixed sequence
    def debiobrick(self):
# the function first checks if the sequence is a BioBrick with an if clause
        if self.sequence[:20] == "GAATTCGCGGCCGCTTCTAG":
            if self.sequence[:22] == "GAATTCGCGGCCGCTTCTAGAG":
                dbb = self.sequence[22:-21]
            else:
                dbb = self.sequence[20:-21]
            return dbb
        else: print "Not a BioBrick"

# this block parses the arguments and puts the file into a readable format
args = parser.parse_args()
filename = args.filename
seq = SeqIO.read(filename,"fasta")
seq = seq.seq
sequence=Sequence(seq)

# this detects the -b (BioBrick) option being selected and prints the
# BioBrick function
if args.b:
    print sequence.biobrick()
# this detects the -d (deBioBrick) option being selected and prints the
# deBioBrick function
if args.d:
    print sequence.debiobrick()


