#!/usr/bin/python

__author__ = "Peter Williams"
__copyright__ = "Copyright 2016"
__credits__ = ["Peter Williams"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Peter Williams"
__email__ = "peter.williams@yale.edu"

### Usage:      python3 distcalc.py -i <input file> -a <index of residue 1> -b <index of residue 2> -f <output format> -o <.txt output filename>
### Examples:   python3 distcalc.py -i sample-input.pdb -a 3 -b 20 -f 2 -o sample-output.txt 
###				python3 distcalc.py -i sample-input.pdb -a 3 -b 20 -f 1 -o sample-output1.txt  
###				python3 distcalc.py -i sample-input.pdb -a 3 -b 20 -f 1
###				{python3 
###				>>>from distcalc import distance
###				>>>dist=distance('sample-input.pdb',3,20,0,'')
###				} 
### Note:       	Calculates the distance between two C-alpha molecules in the backbone of a polypeptide/protein
###						both -f and -o flags are optional
### Input Format:	
###		-i:	input file:		pdb file of a protein
###		-a:	index1:			integer index of amino acid in sequence to calculate the distance of its alpha C to anothers
###		-b:	index2:			integer index of other amino acid in sequence to calculate the distance of its alpha C to that of index1
###		-f:	outopts: 		integer (either 0,1,2) to determine output format, defaults to outputting in line
###							0 : outputs distance in angstroms to command line
###							1 : outputs distance in angstroms to txt output file
###							2 : outputs distance in angstroms to txt output file and amino acid sequence of the protein in the pdb file
###		-o:	outpoutfile:	name of txt file for output options 1 and 2, defaults to "output.txt"

### Import libraries
import argparse
import numpy as np

### Allocate input arguments
if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='Euclidian Distance Calculator between Alpha Carbons in pdb file')
	parser.add_argument('-i', '--inputfile', help='input file name', required=True)
	parser.add_argument('-a', '--index1', help='index of first residue in sequence', required=True)
	parser.add_argument('-b', '--index2', help='index of second residue in sequence', required=True)
	parser.add_argument('-f', '--outopts', nargs='?',default=0, help='value corresponding to desired output')
	parser.add_argument('-o', '--outputfile', nargs='?',default='output.txt', help='output file name')
	args = parser.parse_args()
	#distance(args.inputfile, args.index1, args.index2,args.outopts, args.outputfile)

### Implementation
def distance(filename, index1, index2, outopts, outname):
	#Convert non-string inputs into their proper form
	index1=int(index1)
	index2=int(index2)
	outopts=int(outopts)

	######################################################################################
	### Initialize variables 
	f=open(filename,'r')
	aacount=0 #index of current residue/amino acid as run through pdb file
	dRvec=np.zeros(3) # difference vector between alpha carbons of interest
	addedvec=0 # index that keeps track how many of the alpha carbons have 
				 #been accounted for in dRvec (when =2 have both alpha carbons)
	resname='' # string containing the 3 letter abreviation of the residue
	if outopts==2:
		resnames='' # if keeping track of AA sequence, introduce string resnames to log them

	########################################################################################
	### Read through pdb file and log information
	# Run through line by line in pdb file looking for coordinates of interest and then create
	#     difference vector dRvec
	for x in f:
		# If the line in the pdb file contains coordinates it starts with Atom and has the following form
		# x[0:4] ="Atom"
		# x[6:11] = Atom serial number
		# x[12:16] = Atom Name
		# x[16] = Alternate location indicator
		# x[17:20] = Residue Name
		# x[21] = Chain Identifier
		# x[22:26] = Residue Sequence number
		# x[26] = Code for insertions of residues
		# x[30:38] = x orthogonal angstrom coordinate
		# x[38:46] = y orthogonal angstrom coordinate
		# x[46:54] = z orthogonal angstrom coordinate
		
		# Check if at portion of pdb file that lists coordinates
		if x[0:4]=="ATOM":
			# Initialize a count of the index of the amino acid currently sorting through
			#    Ignore atoms that are not part of amino acid residues 
			#		(ignore standard nucleic acids as they do not have Alpha Cs)
			if aacount==0 and x[17]!=' ':
				aacount+=1
				resname=x[17:20];
				if outopts==2:
					resnames+=resname+', ' 
			# Check if at a new residue. If you are, increase the count of amino acids in sequence,
			#    update resname to reflect that of the current amino acid, and add resname to resnames
			#    if it is desired to print out the amino acid sequence
			if x[17:20]!=resname and x[17]!=' ':
				aacount+=1
				resname=x[17:20]
				if outopts==2:
					resnames+=resname+', '
			# Check if the residue is the residue of the desired indices
			#    Look for alpha carbon and log coordinates into dRvec
			#    Account for the logging of coordinates with addedvec
			#    Choose convention: dR=R[index1]-R[index2]
			#    Distance between is |dR|=np.linalg.norm(dRvec)
			if aacount==index1:
				if x[12:16]==" CA ":
					dRvec[0]+=float(x[30:38])
					dRvec[1]+=float(x[38:46])
					dRvec[2]+=float(x[46:54])
					addedvec+=1
			if aacount==index2:
				if x[12:16]==" CA ":
					dRvec[0]+= -float(x[30:38])
					dRvec[1]+= -float(x[38:46])
					dRvec[2]+= -float(x[46:54])
					addedvec+=1
		# Exit for loop if have difference vector and do not desire to 
		#    print the entire AA sequence
		if outopts!=2 and addedvec==2:
			break
	# Done reading pdb file, close it out
	f.close()

	###########################################################################################################
	### Output
	# Format output in the desired form
	if outopts==0:
		return np.linalg.norm(dRvec)
	else:
		fout=open(outname,'w')
		if outopts==1:
			fout.write('Distance between Alpha Carbon on residue %d from Alpha Carbon on residue %d = %f Angstroms'%(index1, index2, np.linalg.norm(dRvec)))
		if outopts==2:
			fout.write('Distance between Alpha Carbon on residue %d from Alpha Carbon on residue %d = %f Angstroms'%(index1, index2, np.linalg.norm(dRvec)))
			fout.write('\n')
			fout.write('\n')
			fout.write('Amino Acid Sequence: \n')
			fout.write(resnames[0:-2])
		fout.close()
	
### Run
if __name__ == '__main__':
	distance(args.inputfile, args.index1, args.index2,args.outopts, args.outputfile)

