#@file: process_kikn_info.py
#@brief: Process kink data from KinkFinder
#@author: Rebecca F. Alford

import sys, os
import argparse
import numpy as np
from os import listdir
from os.path import isfile, join
import csv

# Headers for kink file csvs
# pdb_code, Helix_Start, Helix_End, Kink_Position,Kink_Start,Kink_End,Kink_Angle,sequence,n_radius,n_rmsd,c_radius,c_rmsd,Estimated_Error

path = "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/helix-kinks/kink_files/3f7v"
kink_files = [ f for f in listdir(path) if isfile(join(path, f)) ]

data_for_pdb = []
for kink_f in kink_files: 

	model_number = kink_f.split(".")[0].split("_")[1]
	with open( kink_f ) as csvfile: 
			csvreader = csv.DictReader(csvfile, delimiter=',')

			helix_start = 0
			helix_end = 0
			helix_num = 1

			already_observed = []

			for row in csvreader: 

				# Collect data for this helix
				pdb_code = row['pdb_code'][0:4]
				curr_helix_start = int(row['Helix_Start'])
				curr_helix_end = int(row['Helix_End'])
				curr_kink_position = row['Kink_Position']
				curr_kink_angle = row['Kink_Angle']

				if ( (helix_start != curr_helix_start) and (helix_end != curr_helix_end) ):
					
					if ( not (curr_helix_start in already_observed) ): 

						# store the data as a string
						datastr = pdb_code + " " + model_number + " " + str(helix_num) + " " + str(curr_helix_start) + " " + str(curr_helix_end) + " " + curr_kink_position + " " + curr_kink_angle
						data_for_pdb.append( datastr )
						helix_start = curr_helix_start
						helix_end = curr_helix_end
						helix_num += 1
						already_observed.append( helix_start )


# Write the data to an output file
with open( "test_output.csv", 'w' ) as f: 
	f.write( "pdb_code model_number helix_num helix_start helix_end kink_position kink_angle\n" )
	for line in data_for_pdb: 
		f.write( line + "\n" )

