# assemblywizard
##Synopsis
Assembly wizard is a suite of programs designed to automate the construction of DNA sequences.

The 3 programs are:
(1) PySlice
(2) biobricker
(3) paperclipper

They are all designed to be run in the linux terminal.
##Code Example
###PySlice

PySlice is designed to import a sequence as a fasta file and run a restriction digest on it with either a pre-defined or user defined restriction enzyme. This is useful for looking at the sticky ends produced in each fragment. It also tells you the strand length so DNA examples can be identified using gel electrophoresis.

The program is run in the terminal by typing python before it.

The command should be typed as:

```
python PySlice.py <filename>, <enzyme> [-d,-e].
```

The file must be a fasta file.

The enzyme should be typed in this space e.g. "EcoRI" or "NotI" (without quotes). If you would like to set your own enzyme type "new".

If you type "new" you will be asked for the name of the enzyme, the recognition pattern which should be in the order on the sense strand and the cut index. The cut index is defined as the base immediately before the cut on the sense strand if the recognition site is numbered left to right, 1 being the first base in the recognition site.


e.g. EcoRI, G  | A A T T C - the cut index is 1

e.g. PstI, C T G C A | G - the cut index is 5

The options in the command line are:
`-d`  Digest sequence
`-r`  Prints the locations of recognition sites.

Example:

An example use of PySlice is made using the "RestrictionExample.fasta" file which is a fasta file that contains random DNA bases. If you wanted to run a digest on the sequence using EcoRI you would type:
```
python PySlice.py RestrictionExample.fasta EcoRI -d
```

If you wanted to know just where the restriction sites were you would type:
```
python PySlice.py RestrictionExample.fasta EcoRI -r
```

To use your own restriction enzyme type:
```
python PySlice.py RestrictionExample.fasta new -d
```

The terminal will prompt you to input the information for a new enzyme.

New enzymes can be manually added to 'renz' dictionary on line 141 in the format `"key":"enzymename","recognitionsequence",cutindex`. This is useful if you are using an enzyme repeatedly.

###biobricker

biobricker is a script for producing biobricks from DNA sequences and converting biobricks to normal DNA sequences. It is based on the RFC10 recommendation from the BioBrick Foundation for making universally compatible DNA parts. More information can be found at: [http://parts.igem.org/Help:Standards/Assembly/RFC10]

The program is run in the linux terminal by writing python before it.

The command should be typed as:

```
python biobricker.py <filename> [-b,-d]
```

The file must be a fasta file.

The options in the command line are:

`-b` Produces a BioBrick from a DNA sequence
`-d` Produces a normal DNA sequence from a BioBrick

Example:

To convert a DNA sequence into a BioBrick we will use the "rhodopsin.fasta" sequence which we want to convert to a BioBrick. To do this we type:

```
python biobricker.py rhodopsin.fasta -b
```

This outputs the sequence as a BioBrick into the terminal.

To convert the other way from BioBrick to normal sequence we will use the "biobrick.fasta" file which contains the rhodopsin sequnce with a BioBrick prefix and suffix. To do this type:
```
python biobricker.py biobrick.fasta -d
```

This ouputs the sequence without the prefix and suffix.

###paperclipper

paperclipper is a script for working out the oligonucleotides needed for performing a PaperClip assembly. PaperClip assembly is an assembly method designed by Trubitsyna et al. for the rapid assembly of DNA libraries. Oligonucleotides are required for each part to make clips which join the parts together in any order. More information about the process is avaiable at: [https://doi.org/10.1093/nar/gku829]. This tool is designed to speed up the process of making the olignucleotides.

The program is run in the linux terminal by writing python before it:

The command should be typed as:

```
python paperclipper.py <filename>
```

The file must be a fasta file.

This outputs the sequences for the oligonucleotides.

Example:

To produce the oligonucleotides for a sequence we will use the "rhodopsin.fasta" file. To do this type:
```
python paperclipper.py rhodopsin.fasta
```

The output prints the oligonucleotides that need to be produced. It also shows what the half clips will look like once the four olignucleotides are annealed.

##Installation
The scripts can be found at [https://github.com/alihum/assemblywizard].


