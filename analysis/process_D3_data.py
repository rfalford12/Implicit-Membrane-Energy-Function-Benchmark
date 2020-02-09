#!/usr/bin/env python
#@file: process_docking_data.py
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Parse files containing interface RMS values and dG values

import sys, os, subprocess
import argparse
import numpy as np

def remove_extra_data( filename, exclude  ): 

	# Preserve the header
	raw_header = subprocess.check_output( ['head', '-n2', filename ] )
	header = raw_header.split("\n")[1]

	# Remove extra headers
	os.system( "grep -v '^SEQUENCE' " + filename + " > temp.sc")
	os.system( "grep -v '^" + header + "' temp.sc > temp2.sc" )

	# Add the header back to the file
	with open( "headerfile", 'w' ) as f: 
		f.write( header + "\n" )

	os.system( "cat headerfile temp2.sc > temp3.sc" )

	# Read in the resulting file
	with open( "temp3.sc", 'rt' ) as g: 
		contents = g.readlines()
		contents = [ x.strip() for x in contents ]

	# Iterate through lines in the file, make a list of excluded models
	included = []
	exclusion_list = []
	for x in range(0, len(contents)): 

		contents[x] = contents[x].replace( "_0001", "" )
		dataline = contents[x].split(" ")
		included.append(contents[x])

	# Clean up after yourself
	os.system( "rm temp.sc" )
	os.system( "rm temp2.sc" )
	os.system( "rm temp3.sc" )
	os.system( "rm headerfile" )

	return included, exclusion_list

def get_duplicates_list( docking_models, angles_for_models ):

	# Get list of descriptions
	docking_model_dsc = []
	angles_for_model_dsc = []
	for model in docking_models: 
		docking_model_dsc.append( model.split(" ")[-1] )
	for model in angles_for_models: 
		angles_for_model_dsc.append( model.split(" ")[-1] )

	# Get a list of descriptions that are not included in both sets
	diff_one_way = list(set(docking_model_dsc) - set(angles_for_model_dsc))
	diff_other_way = list(set(angles_for_model_dsc) - set(docking_model_dsc))
	full_diff = diff_one_way + diff_other_way 
	return full_diff

def remove_exclusive_and_duplicate_lines( model_list, exclusion_list ): 

	included_models = []
	seen_descripttion = []
	for model in model_list: 

		description = model.split(" ")[-1]
		if ( description in seen_descripttion ): 
			print "Duplicate description, " + description + " skipping..."
		else: 
			seen_descripttion.append( description )
			if ( description in exclusion_list ): 
				print "Excluding " + description + "..."
			else: 
				included_models.append( model )

	return included_models

def sort_models_by_description( filename ) : 

	# Read in each line, split by the delimiter
	with open( filename, 'rt' ) as f: 
		contents = f.readlines()
		contents = [ x.strip() for x in contents ]
		scoredata = [ x.split(" ") for x in contents ]

	# Grab the descriptions
	temp = [ scoredata[x][-1].split("_")[-1] for x in range(0, len(scoredata)) ]
	description_numbers = []
	for desc in temp: 
		if ( desc.isdigit() ): 
			description_numbers.append( int(desc) )
	args = np.argsort(description_numbers)

	# Make the name of a sorted file
	sorted_file = filename.split(".")[0] + "_sorted.sc"
	with open( sorted_file, 'wt' ) as f: 
		for arg in args: 
			f.write( contents[arg] + "\n" )
	f.close()

	return sorted_file



def write_models_to_file( model_list, filename ): 

	# Write out the parsed PDB file
	with open( filename, 'w' ) as f: 
		for x in model_list: 
			f.write( x + "\n" )

def main( args ): 

	parser = argparse.ArgumentParser(description="Return a parsed, combined file with interface and angle data.")
	parser.add_argument("--basename", type=str, help="Basename for output files")
	parser.add_argument("--local", type=bool, help="Local docking option")
	input_args = parser.parse_args()

	# Get file basename
	basename = input_args.basename

	# Store input filenames
	docked_models_file = basename + "_docked.sc"
	interfaces_file = basename + "_interfaces.sc"
	if ( input_args.local ): 
		docked_models_file = basename + "_local_refine.sc" 
		interfaces_file = basename + "_interfaces_local.sc"

	# Remove extra data from angles file
	docked_models_data, exclusion_list = remove_extra_data( docked_models_file, True )
	interface_data, exclusion_list2 = remove_extra_data( interfaces_file, True )

	# Get the complete exclusion list, accounting for duplicates
	no_duplicates = get_duplicates_list( docked_models_data, interface_data )
	full_excluded = no_duplicates + exclusion_list

	# Get pruned list of docking models
	parsed_docking_models = remove_exclusive_and_duplicate_lines( docked_models_data, full_excluded )
	parsed_interfaces_for_models = remove_exclusive_and_duplicate_lines( interface_data, full_excluded )

	# Write the resulting models to file
	parsed_models = basename + "_docked_parsed.sc"
	parsed_interfaces = basename + "_interfaces_parsed.sc"
	write_models_to_file( parsed_docking_models, parsed_models )
	write_models_to_file( parsed_interfaces_for_models, parsed_interfaces )

	# Before sorting, grap the first few columns to make life easier
	os.system( "awk '{print $2,$3,$4,$5,$6,$7,$8,$NF}' " + parsed_models + " > parsed_models.temp")

	# Sort the parsed files
	sorted_parsed_interface = sort_models_by_description( parsed_interfaces )
	sorted_parsed_models = sort_models_by_description( "parsed_models.temp" )

	# Post-process and paste the files together
	os.system( "awk '{print $6,$7,$8,$9,$10}' " + parsed_interfaces + " > parsed_int.temp" )
	if ( input_args.local ): 
		os.system( "paste parsed_models.temp parsed_int.temp > " + basename + "_data_local.sc" )
	else: 
		os.system( "paste parsed_models.temp parsed_int.temp > " + basename + "_data.sc" )
	os.system( "rm parsed_int.temp" )
	os.system( "rm parsed_models.temp" )
	os.system( "rm " + parsed_interfaces )
	os.system( "rm " + parsed_models)
	os.system( "rm " + sorted_parsed_models )



if __name__ == "__main__" : main(sys.argv)



