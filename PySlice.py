# -*- coding: utf-8 -*-
"""
Created on Sat Mar  4 13:24:11 2017

@author: ali
"""
# Imports necessary modules "sys" for dealing with command line arguments
# "SeqIO" for reading a fasta sequence file
import argparse
from Bio import SeqIO

# makes a restriction enzyme class
class RestrictionEnzyme:
# the __init__ initialises the instance of a RestrictionEnzyme object
#  which takes the variable of the name, the sequence of the recognition 
# site (cutsite), and the cutindex which is the number of the base immediately
# before the cut numbered from the start of the recognition sequence
    def __init__(self,name,cutsite,cutindex):
        self.name = name
        self.cutsite = cutsite
        self.cutindex = cutindex
        self.reverse_cutsite = cutsite[::-1]
        self.cutsite_length = len(cutsite)

# a printable string that tells you the restriction enzme and its recognition
# sequence     
    def __str__(self):
        return "This {} restriction enzyme with recognition sequence {} ."\
        .format(self.name,self.cutsite)
    
# a function for producing the complementary sequence to the sequence variable
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
    
# calculates the points at which the sense strand will be cut and compiles
# them into a list
    def cutpoints(self,sequence):
        x = []
        for i in range(0, len(sequence)-len(self.cutsite)):
            if sequence[i:i+len(self.cutsite)] == self.cutsite:
                x.append(i+self.cutindex)
        return x
   
# this function finds the locations of each recognition sequence and returns 
   #it as as list
    def findrecognitionsites(self,sequence):
        cp = self.cutpoints(sequence)
        return cp
    
# this function prints the output of findrecognitionsites function
    # if the list is emtpy it prints a message
    def printrecognitionsites(self,sequence):
        cp2 = self.findrecognitionsites(sequence)
        if cp2 == []:
            print("No recognition site found.")
        else:
            for i in cp2:
                    print "{} restriction site from {} to {} bps"\
                .format(self.name,i,i+len(self.cutsite)-1)
    
# a functions that digests the sequence with a specific restriction enzyme
# and prints how many fragments, the length  of each fragment and a graphical
# display of the sticky ends    
    def restrictionfragments(self,sequence):
         #makes a list of the points where the sequence is to be cut
         rs = self.cutpoints(sequence)
         # this is the length of the overhang calculated from the 
         cutlength = (len(self.cutsite)-self.cutindex*2)
         rs.append(len(sequence))
         complement = self.complement(sequence)
         print "There will be {} fragment(s)".format(len(rs))
         count = 1
         
         if cutlength >= 0:
             for i in range(len(rs)):
                 # count makes sure the first fragment is displayed properly
                 if count == 1:
                     print "Fragment 1 ({}bps):\n".format(rs[i])
                     print  sequence[:5] + "."*5 + sequence[rs[i]-5:rs[i]]
                     print "|"*5 + " "*5 + "|"*5
                     print complement[:5] + "."*5 + \
                     complement[rs[i]-5:rs[i]+cutlength]                 
                     count = count + 1
                 else:
                     print "\nFragment {}({}bps):\n"\
                     .format(count,rs[i]-rs[i-1])
                     print sequence[rs[i-1]:rs[i-1]+cutlength+5] + \
                     "."*len(self.cutsite) +sequence[rs[i]-5:rs[i]]
                     print " "*(cutlength) + "|"*5 + " "*len(self.cutsite) +\
                     "|"*5
                     print " "*(cutlength) +complement[rs[i-1]+\
                     cutlength:rs[i-1]+cutlength+5] + "."*len(self.cutsite)\
                     + complement[rs[i]-5:rs[i]+cutlength]
                     count = count + 1
         else:
             #if the cutindex is more than half the recognition site this else
         #clause makes sure that it is displayed properly graphically
            cutlength = abs(cutlength)
            for i in range(len(rs)):
                 if count == 1:
                     print "Fragment 1 ({}bps):\n".format(rs[i])
                     print  sequence[:5] + "."*5 + sequence[rs[i]-5:rs[i]\
                     +cutlength]
                     print "|"*5 + " "*5 + "|"*5
                     print complement[:5] + "."*5 + \
                     complement[rs[i]-5:rs[i]]                 
                     count = count + 1
                 else:
                     print "\nFragment {}({}bps):\n"\
                     .format(count,rs[i]-rs[i-1])
                     print " "*(cutlength) +sequence[rs[i-1]:rs[i-1]+5] \
                     + "."*len(self.cutsite) +sequence[rs[i]-5-cutlength:rs[i]]
                     print " "*(cutlength) + "|"*5 + " "*len(self.cutsite)\
                     + "|"*5
                     print complement[rs[i-1]-cutlength:rs[i-1]+5] \
                     + "."*len(self.cutsite) \
                     + complement[rs[i]-5-cutlength:rs[i]-cutlength]
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
renz = {
"EcoRI":("EcoRI","GAATTC",1),
"NotI":("NotI","GCGGCCGC",2),
"NheI":("NheI","GCTAGC",1),
"PstI":("PstI","CTGCAG",5),
"SpeI":("SpeI","ACTAGT",1),
"EcoRV":("EcoRV","GATATC",3),
"XbaI":("XbaI","TCTAGA",1),
}

parser = argparse.ArgumentParser(description="PySlice is a program for "\
+ "manipulating nucelotide sequences with restriciton enzymes")

parser.add_argument('filename',
    help="Insert name of fasta file - must be in same directory" )
    
parser.add_argument('restrictionenzyme',choices=renz.keys(),
    help="Choose the restriction enzyme to use")
    
parser.add_argument('-d',action='store_true')

parser.add_argument('-r',action='store_true')

args = parser.parse_args()
file = args.filename
seq = SeqIO.read(file,"fasta")
seq = seq.seq

re = args.restrictionenzyme
x = RestrictionEnzyme(renz[re][0],renz[re][1],renz[re][2])   

if args.d:
    print x.restrictionfragments(seq)
    
if args.r:
    print x.printrecognitionsites(seq)

            
