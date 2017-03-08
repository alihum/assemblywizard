# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:24:11 2017

@author: ali
"""
import sys
from Bio import SeqIO

file = sys.argv[1]
seq = SeqIO.read(file,"fasta")
seq = seq.seq

class RestrictionEnzyme:
    
    def __init__(self,name,cutsite,cutindex):
        self.name = name
        self.cutsite = cutsite
        self.cutindex = cutindex
        self.reverse_cutsite = cutsite[::-1]
        self.cutsite_length = len(cutsite)
    
    def __str__(self):
        return "This is the {} restriction enzyme with cutsite {}"\
        .format(self.name,self.cutsite)
    
    def complement(self,sequence):
        complement = ""
        for i in range(len(sequence)):
            if sequence[i] == "A":
                complement = complement + "T"
            elif sequence[i] == "T":
                complement = complement + "A"
            elif sequence[i] == "C":
                complement = complement + "G"
            else:complement = complement + "C"
        return(complement)
    
    def cutpoints(self,sequence):
        x = []
        for i in range(0, len(sequence)-len(self.cutsite)):
            if sequence[i:i+len(self.cutsite)] == self.cutsite:
                x.append(i+self.cutindex)
        return x
    
    def findrecognitionsites(self,sequence):
        cp = self.cutpoints(sequence)
        return cp
    
    def printrecognitionsites(self,sequence):
        cp2 = self.findrecognitionsites(sequence)
        if cp2 == []:
            print("No recognition site found.")
        else:
            for i in cp2:
                    print "{} restriction site from {} to {} bps"\
                .format(self.name,i,i+len(self.cutsite)-1)
    
    def restrictionfragments(self,sequence):
         rs = self.cutpoints(sequence)
         cutlength = (len(self.cutsite)-self.cutindex)-1
         rs.append(len(sequence))
         complement = self.complement(sequence)
         print "There will be {} fragment(s)".format(len(rs))
         count = 1
         for i in range(len(rs)):
             if count == 1:
                 print "Fragment 1 ({}bps):".format(rs[i])
                 print  sequence[:3] + "."*3 + sequence[rs[i]-3:rs[i]]
                 print "|"*3 + " "*3 + "|"*3
                 print complement[:3] + "."*3 + \
                 complement[rs[i]-3:rs[i]+cutlength]                 
                 count = count + 1
             else:
                 print "Fragment {}({}bps):"\
                 .format(count,rs[i]-rs[i-1])
                 print sequence[rs[i-1]:rs[i-1]+\
                 len(self.cutsite)-self.cutindex+3] + "."*len(self.cutsite) + \
                 sequence[rs[i]-3:rs[i]]
                 print " "*(cutlength) + "|"*4 + " "*len(self.cutsite) + "|"*3
                 print " "*(cutlength) + \
                 complement[rs[i-1]+cutlength:rs[i-1]+cutlength+4] \
                 + "."*len(self.cutsite) + complement[rs[i]-3:rs[i]+cutlength]
                 count = count + 1
    
    # prints the overhang that is produced from a restriction enzyme
    def overhang(self):
        print self.cutsite[:self.cutindex] + "      " \
        + self.cutsite[self.cutindex:]
        print self.reverse_cutsite[:self.cutsite_length-self.cutindex] \
        + "      " + self.reverse_cutsite[self.cutindex]

# Place to store restriction enzymes
# They are inserted in the format (Name,Recognition Sequence,CutIndex)
# The cut index is the number of the base immediately before where
# restriction occurs e.g. G AATTC is an index of 1 | GC GGCCGC = 2
# User defined restriction enzymes can be stored here
EcoRI = RestrictionEnzyme("EcoRI","GAATTC",1)
NotI = RestrictionEnzyme("NotI","GCGGCCGC",2)
NheI = RestrictionEnzyme("NheI","GCTAGC",1)
PstI = RestrictionEnzyme("PstI","CTGCAG",5)
SpeI = RestrictionEnzyme("SpeI","ACTAGT",1)
EcoRV = RestrictionEnzyme("EcoRV","GATATC",3)
XbaI = RestrictionEnzyme("XbaI","TCTAGA",1)

# Detects the "digest" options in the command line and executes this
if sys.argv[2]=="digest":
# Gives user the option to input their own restriction enzyme or use of the
# built in ones    
    if raw_input("Use own restriction enzyme? Y/N : ") == "Y":
        x = raw_input("Type the name of restriction enzyme: ")
        y = raw_input("Type recognition sequence: ")
        z = int(raw_input("Type cut index: "))
        re = RestrictionEnzyme(x,y,z)
        print re.restrictionfragments(seq)
    else:
        x = raw_input("Type the name of a restriction enzyme (e.g. EcoRI): ")
        if x == "EcoRI":
            print EcoRI.restrictionfragments(seq)
        elif x == "NotI":
            print NotI.restrictionfragments(seq)
        elif x == "NheI":
            print NheI.restrictionfragments(seq)
        elif x == "PstI":
            print PstI.restrictionfragments(seq)
        elif x == "SpeI":
            print SpeI.restrictionfragments(seq)
        elif x == "EcoRV":
            print EcoRV.restrictionfragments(seq)
        elif x == "XbaI":
            print XbaI.restrictionfragments(seq)
            
            