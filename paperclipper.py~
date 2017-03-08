# -*- coding: utf-8 -*-
"""
Created on Sun Mar  5 15:58:20 2017

@author: ali
"""

import sys
from Bio import SeqIO
#
fasta_file = sys.argv[1]

seq = SeqIO.read(fasta_file,"fasta")
#
seq = SeqIO.read("id120009419.seq","fasta")

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