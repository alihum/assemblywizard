#!/usr/bin/python
"""
PySlice is designed to import a sequence as a fasta file and run a restriction\
 digest on it with either a pre-defined or user defined restriction enzyme. \
 The program only works with Type II restriction enzymes. This is useful for \
 looking at the sticky ends produced in each fragment. It also tells you the \
 strand length so DNA examples can be identified using gel electrophoresis.

The program is run in the terminal by typing python before it.

The command should be typed as:

python PySlice.py <filename>, <enzyme> [-d,-e].

The file must be a fasta file.

The enzyme should be typed in this space e.g. "EcoRI" or "NotI"\
 (without quotes). If you would like to set your own enzyme type "new".

If you type "new" you will be asked for the name of the enzyme, the \
recognition pattern which should be in the order on the sense strand and the \
cut index. The cut index is defined as the base immediately before the cut on \
the sense strand if the recognition site is numbered left to right, 1 being\
 the first base in the recognition site.


e.g. EcoRI, G  | A A T T C - the cut index is 1

e.g. PstI, C T G C A | G - the cut index is 5

The options in the command line are:
`-d`  Digest sequence
`-r`  Prints the locations of recognition sites.
"""
# Imports necessary modules "sys" for dealing with command line arguments
# "SeqIO" for reading a fasta sequence file
import argparse
from Bio import SeqIO

# makes a restriction enzyme class
class RestrictionEnzyme:
    """Form a restriction enzyme class:
    the __init__ initialises the instance of a RestrictionEnzyme object
    which takes the variable of the name, the sequence of the recognition 
    site (cutsite), and the cutindex which is the number of the base \
    immediately
    before the cut numbered from the start of the recognition sequence
    """  
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
                     print "Fragment 1 ({}bps):".format(rs[i])
                     print "0bp to {}bp\n".format(rs[i])
                     print  sequence[:5] + "."*5 + sequence[rs[i]-5:rs[i]]
                     print "|"*5 + " "*5 + "|"*5
                     print complement[:5] + "."*5 + \
                     complement[rs[i]-5:rs[i]+cutlength]                 
                     count = count + 1
                 else:
                     print "\nFragment {}({}bps):"\
                     .format(count,rs[i]-rs[i-1])
                     print "{}bp to {}bp\n".format(rs[i-1],rs[i])
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
                     print "Fragment 1 ({}bps):".format(rs[i])
                     print "0bp to {}bp\n".format(rs[i])
                     print  sequence[:5] + "."*5 + sequence[rs[i]-5:rs[i]\
                     +cutlength]
                     print "|"*5 + " "*5 + "|"*5
                     print complement[:5] + "."*5 + \
                     complement[rs[i]-5:rs[i]]                 
                     count = count + 1
                 else:
                     print "\nFragment {}({}bps):"\
                     .format(count,rs[i]-rs[i-1])
                     print "{}bp to {}bp\n".format(rs[i-1],rs[i])
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

# Dictionary to store restriction enzymes
# They are inserted in the format (Name,Recognition Sequence,CutIndex)
# The cut index is the number of the base immediately before where
# restriction occurs e.g. G AATTC has an index of 1 | GC GGCCGC = 2
# User defined restriction enzymes can be stored here
renz = {
"new":(),
"EcoRI":("EcoRI","GAATTC",1),
"NotI":("NotI","GCGGCCGC",2),
"NheI":("NheI","GCTAGC",1),
"PstI":("PstI","CTGCAG",5),
"SpeI":("SpeI","ACTAGT",1),
"EcoRV":("EcoRV","GATATC",3),
"XbaI":("XbaI","TCTAGA",1),
}

#Command-line argument parser so the program can be run from the linux shell
parser = argparse.ArgumentParser(description="PySlice is designed to import\
 sequence as a fasta file and run a restriction digest on it with either a \
 pre-defined or user defined restriction enzyme. The program only works with \
 Type II restriction enzymes. This is useful for looking at the sticky ends\
 produced in each fragment. It also tells you the strand length so DNA\
 examples can be identified using gel electrophoresis.")

#The first two arguments are positional arguments which take the file and the
# restriction enzyme to be used
parser.add_argument('filename',
    help="Insert name of fasta file - must be in same directory" )
#The restriction enzyme can be one of the preset ones or "new" which allows 
    #user to define their own enzyme
parser.add_argument('enzyme',choices=renz.keys(),
    help="Choose the restriction enzyme to use")    
#This adds an option to run a digest on the sequence
parser.add_argument('-d',action='store_true',
    help="Run a digest on sequence")
    #This adds an option to print the locations of the recognition sequences
parser.add_argument('-r',action='store_true',
    help="Print the locations of restriction fragments")

#This piece of code parses the file using SeqIO also has an error handler if
# file is not recognised
try:
    args = parser.parse_args()
    file = args.filename
    seq = SeqIO.read(file,"fasta")
    seq = seq.seq
except IOError: print "No such file!"
#The if clause checks to see if the user wants to use a pre-defined enzyme
# or a new enzyme by looking for "new" in the enzyme argument
if args.enzyme == "new":
    #set of raw_input statements asks for user input to create an instance
    # the RestrictionEnzyme class based on their input
    name = raw_input("Type enzyme name (e.g. 'EcoRI'):")
    recognition = raw_input("Type recognition sequence (e.g. 'GAATTC'):")
    cutindex = int(raw_input("Type cut index (e.g. 1):"))
    x = RestrictionEnzyme(name,recognition,cutindex)
else:
    #if a pre-defined enzyme is desired this else clause creates an instance
    #of the RestrictionEnzyme class based on the renz dictionary which
    #stores pre-defined restriction enzymes    
    re = args.enzyme
    x = RestrictionEnzyme(renz[re][0],renz[re][1],renz[re][2])   

# when the -d option is selected this prints the restriction fragments
# function
if args.d:
    print x.restrictionfragments(seq)

# when the -r option is selected this prints the recognition sites function        
if args.r:
    print x.printrecognitionsites(seq)
