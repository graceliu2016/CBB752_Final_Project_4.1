**README for CBB_Bioinformatics_FinalProject_4.1**
---------------------------------------------------------------
Tool that calculates the distance between two alpha carbons from a PDB file

# The python tool that accomplishes this task is named distcalc.py
distcalc.py takes 3 required inputs (inputfile, index1, and index2) and has two optional (outopts, and outputfile)
The indices (index1 and index2) correspond to the alpha cabon of the nth residue. The input inputfile is the name of a corresponding pdb file from which to calculate a distance.
The tool can be called from both a terminal or inside python with the formatting listed below.

Usage:      python3 distcalc.py -i < input file > -a < index of residue 1 > -b < index of residue 2 > -f < output format > -o <.txt output filename>

Examples:
```{r NCBI_python, engine="python", highlight=TRUE}
# Usage from terminal:
	
    python3 distcalc.py -i sample-input.pdb -a 3 -b 20 -f 2 -o testout.txt 
     
    python3 distcalc.py -i sample-input.pdb -a 3 -b 20 -f 1 -o testout.txt
        
    python3 distcalc.py -i sample-input.pdb -a 3 -b 20 -f 1
            
# In line usage in python:
  	
    python3 
       		
	>>>from distcalc import distance
	
  	>>>dist=distance('1a3n.pdb',3,20,0,'')
```

Input Formats:
```{r NCBI_python, engine="python", highlight=TRUE}
	-i:	input file:		pdb file of a protein
	-a:	index1:			integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
	-b:	index2:			integer index of other amino acid in sequence to calculate the distance of its alpha C to that of index1
	-f:	outopts: 		integer (either 0,1,2) to determine output format, defaults to outputting in line
						0 : outputs distance in angstroms to command line
						1 : outputs distance in angstroms to txt output file
						2 : outputs distance in angstroms to txt output file and amino acid sequence of the protein in the pdb file
	-o:	outputfile:	name of txt file for output options 1 and 2, defaults to "output.txt"
```

# The R tool that accomplishes this task is named distcalc.rmd

## General

The R piece of code will read a pdb file and find the distance between two selected alpha carbons.

The R tool that accomplishes this task is named distcalc. Distcalc takes 5 required inputs: inputfile, index1, index2, outopts, and outputfile. The indices (index1 and index2) correspond to the alpha carbon of the nth residue. The input inputfile is the name of a corresponding pdb file from which to calculate a distance. The code will run after loading the function from the distcalc.rmd file.

### Input Options for discalc:

1. inputfile - the name of the input file (1a3n.pdb)
2. index1 - integer index of 1st amino acid in sequence to calculate the distance of its alpha C to other (2)
3. index2 - integer index of 2nd amino acid in sequence to calculate the distance of its alpha C to other (9)
4. outops - integer (either 0,1,2) to determine output format
  + 0 : outputs distance in angstroms to console
  + 1 : outputs distance in angstroms to txt output file called output.txt
  + 2 : outputs distance in angstroms to txt output file and amino acid sequence of the protein in the pdb file
5. outputfile - the desired name of the output file ("output.txt")

### Usage

```{r}
distcalc(inputfile = "1a3n.pdb", index1 = 43, index2 = 44, outops = 0, outputfile = "output.txt") 
```
