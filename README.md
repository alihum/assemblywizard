# assemblywizard

Assembly wizard is a suite of programs designed to automate the construction of DNA sequences.

The 3 programs are:
(1) PySlice
(2) biobricker
(3) paperclipper

They are all designed to be run in the linux terminal.

_______
PySlice
_______

PySlice is designed to import a sequence as a fasta file and run a restriction digest on it with either a pre-defined or user defined restriction enzyme.

The program is run in the terminal by typing python before it.

The command should be typed as python PySlice.py <filename>, <enzyme> [-d,-e].

The file must be a fasta file.

The enzyme should be typed in this space e.g. "EcoRI" or "NotI" (without quotes). If you would like to set your own enzyme type "new".

If you type "new" you will be asked for the name of the enzyme, the recognition pattern which should be in the order on the sense strand and the cut index. The cut index is defined as the base immediately before the cut on the sense strand if the recognition site is numbered left to right, 1 being the first base in the recognition site.

            1    2 3 4 5 6
e.g. EcoRI, G  | A A T T C - the cut index is 1

           1 2 3 4 5   6
e.g. PstI, C T G C A | G - the cut index is 5

__________
biobricker
__________


biobricker is a script for producing biobricks from DNA sequences and converting biobricks to normal DNA sequences
