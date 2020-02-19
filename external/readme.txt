What Kink Finder Does
=====================

Any problems, queries, or suggestions, please contact eleanor.law@dtc.ox.ac.uk. 


Kink Finder locates kinks in helices. Based on the local structure of the helix kink, it identifies a residue on the inside of the kink as the kink residue. An angle for each residue is calculated by fitting the coordinates of the 4 backbone atoms of that residue and each of the 5 residues before it (i.e. towards the n-terminus of the protein) to a cylinder, and the backbone atoms of the six residues after it also to a cylinder, and calculating the angle between the two cylinder axes. Kinks are identified starting with the residue with the largest angle, subsequently the residues are examined in decreasing order of angle. A residue is identified as a kink if it is at least 5 residues away from all existing kinks, and there is a residue between it and all existing kinks that has an angle less than 10. 
For each kink, the kink residue is then selected by considering the residue initially identified as the kink, the residue before, and the two after, and choosing the residue most 'inside', based on the 'wobble' angle.
For each angle calculated, an estimation of the error in the angle measurement is calculated based on the quality of fit of the cylinders either side of the kink.

How to install Kink Finder 
============================

Kink Finder requires Python 2.7, which can be obtained from http://www.python.org/
It uses the numpy (http://www.numpy.org/) and argparse (https://pypi.python.org/pypi/argparse) packages.

Linux:
------
Download the file Kink_Finder_linux.tar.gz

tar -zxvf Kink_Finder.tar.gz 

Note:
-----
A pre-compiled C binary for linux is supplied, which performs the cylinder fits. This works on a wide range of linux systems, however it may not work on yours. The source code is supplied, so if it does not work (see troubleshooting section), you may need to compile it on your system, for which instructions are supplied. The previous version of Kink Finder (without error estimation, available at http://www.stats.ox.ac.uk/research/proteins/resources) is provided with binary/executable for mac/windows, but this linux version uses a modified cylinder fitting method.


How to use Kink Finder
======================

For usage instructions type:

python Kink_Finder.py -h


Example Usage:
--------------

python Kink_Finder.py -f ./3EML.pdb -o ./results_eg1/ -d

This will run kink finder on the file 3EML.pdb, and use the helix definitions within the .pdb file. Output files will be written to ./results_eg1/. Summary results for each helix will be displayed to the terminal. This display is reproduced in the helices.csv output file.


Outputs:
--------

helices.csv: Comma separated file, with each line detailing information for an analysed helix: the chain code, the first residue in the helix, the last residue in the helix, the number of the kink residue, and the maximum angle of that kink.

kinks.csv: Comma separated file, with each line reporting any angle over 10 degrees. The 13 residues around the kink residue are shown (6 before, the kink residue itself, and the 6 after).  There are fields for the chain code, the start and end residue number of the helix, the number of the kink residue, the start and end residue number of the 13 residue segment, the kink angle, the sequence of the 13 residues segment around the kink. The next 4 fields give information about the cylinder fits used to calculate the kink angle. These are the rmsd from the cylinder surface, and the radius of the cylinder fitted to the n-terminal and c-terminal side of the kink. The final field is the estimated error which is a function of the sum of the two rmsds.

angles.csv: Comma separated file. Each row starts with the pdb code, chain identifier, and index of the first helix residue concatenated together. Each subsequent entry is the calculated angle for each residue of the helix. The first 5 and last 6 are set to 0, as no angle can be calcuated for these residues.

errors.csv: Comma separated file. Each row starts with the pdb code, chain identifier, and index of the first helix residue concatenated together. Each subsequent entry is the estimated error for the angle at each residue of the helix. The first 5 and last 6 are set to 0, as no angle can be calcuated for these residues.


More advanced usage:
--------------------

To use your own helix definitions, first split the pdb file into its constituent chains, and input the pdb numbers of the first and last residues in each helix, like this:

python Kink_Finder.py -f ./3EMLA.pdb -o ./results_eg2/ -l '7-45 223-270'

For custom helices in multiple chains, edit the helix definitions within the pdb file, which look like this:
HELIX    1   1 SER A    6  ASN A   34  1                                  29

The 5th and 8th entries (which are A in this case) indicate the chain 
The 6th entry is the pdb number of the first residue in the helix
The 8th entry is the pdb number of the last residue in the helix 
The 9th entry indicates the type of helix, in this case 1 meaning alpha helix. Kink Finder will only analyse helices with the code 1 

Be aware that pdb files care about whitespace


Other arguments:
----------------

-d --display display results in terminal whilst Kink Finder is running, with the same information as in the helices.csv file
-p --pymol writes a python pymol script for each helix, named [pdbcode first_residue_of_helix .py] (e.g. 3EMLA7.py), into the directory ./Kink_Finder_pymol_files/ (the name of this directory can be changed in the Kink_Finder.py file). To run this file, cd into a folder containing the coordinate structure, and type:
pymol path_to_python_file.py
Note that cylinders are only rendered as opaque when they are ray traced.  


Further arguments are found in the Kink_Finder.py file, with explanations of their use.


Troubleshooting
===============

Binary/executable errors:
Kink Finder uses a c function to fit its cylinders. A binary (linux) or executable (windows) has been provided. 

"Cylinder binary not in Kink Finder's directory. See readme.txt for further instructions"
The cylinder binary must be in the same directory as the Kink_Finder.py file, otherwise it cannot be called. 

"Cylinder fitting binary not functioning properly. See readme.txt for further instructions"
The pre-compiled binary is not working on your system, so you will need to compile the cylinder binary from the source code provided in the folder ./cylinder_source_code/

To complile (linux):
cd cylinder_source_code/
gcc cylinder.c -o ../cylinder -lm
cd ..

For non-linux systems, c compliers are not installed by default. If your do not have a C Complier currently installed, they can be obtained from a number of places, for example, http://www.mingw.org/ (for windows) and  https://developer.apple.com/xcode/ (for mac), although these are only suggestions, and the author makes no guarentee that these are suitable for your system.

Regardless of your choice of compiler, the -lm flag is the only flag necessary for compilation.

To test the binary, type:

./cylinder 10 0 0 0 0 0 1 -1 0 0 0 1 1 1 0 2 0 -1 3 -1 0 4 0 1 5 1 0 6 0 -1 7 -1 0 8 0 1 9

Which should return the output: 
-0.000000 0.000000 -0.000000 -0.000000 0.000000 1.000000 1.000000 0.000000 0.000000

(The 0.000000 values may be positive or negative)

If these steps do not solve your problems, please email the author (henry.wilman@dtc.ox.ac.uk), and he will endeavour to find a solution for you.


