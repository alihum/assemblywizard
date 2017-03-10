#!/usr/bin/python
"""
paperclipper is a script for working out the oligonucleotides needed for\
 performing a PaperClip assembly. PaperClip assembly is an assembly method \
 designed by Trubitsyna et al. for the rapid assembly of DNA libraries. \
 Oligonucleotides are required for each part to make clips which join the \
 parts together in any order. More information about the process is avaiable \
 at: [https://doi.org/10.1093/nar/gku829]. This tool is designed to speed up\
 the process of making the olignucleotides.

The program is run in the linux terminal by writing python before it:

The command should be typed as:

python paperclipper.py <filename>

The file must be a fasta file.

This outputs the sequences for the oligonucleotides.

Example:

To produce the oligonucleotides for a sequence we will use the\
 "rhodopsin.fasta" file. To do this type:
python paperclipper.py rhodopsin.fasta

The output prints the oligonucleotides that need to be produced.\
 It also shows what the half clips will look like once the four\
 olignucleotides are annealed.

"""

import sys
from Bio import SeqIO

fasta_file = sys.argv[1]

seq = SeqIO.read(fasta_file,"fasta")

# Makes upstream forward oligonucleotide
UF = "GCC" + seq.seq[:40]
# Makes upstream reverse olignucleotide
UR = seq.seq[:37].reverse_complement()

# Selects last 40 bases of part
last_40 = seq.seq[-40:]
#Checks that it does not start with GCC or GCG and continues to remove first
#base until this condition is met
while last_40[:3] == "GCC" or last_40[:3] == "GGC":
    last_40 = seq[1:]
DF = last_40

#Works out downstream reverse olignucelotude which should be 3 bases shorter
# than downstream forward oligo
DR = "GGC" + seq.seq[3-len(DF):].reverse_complement()

# Print statements to graphically display the oligonucleotide sequences
print "\nThese are the oligonucleotide sequences that will make up the half"+ \
"-clips based on the 'Paper-Clip' method of DNA assembly" + \
" (https://doi.org/10.1093/nar/gku829)\n"

print "Upstream forward oligo is:\n\n" + UF
print "\nUpstream reverse oligo is:\n\n" + UR
print "\nDownstream forward oligo is:\n\n" + DF
print "\nDownstream reverse oligo is:\n\n" + DR

print "\nThe upstream half-clip:"
print "\n" + UF
print "   |||||||||||||||||||||||||||||||||||||"
print "   " + UR[::-1]
print "\nThe downstream half-clip:"
print "\n" + DF
print "   |||||||||||||||||||||||||||||||||||||"
print "   " + DR[::-1]
