#!/usr/bin/env python
###################################################################
#@file:         make_refined_models.py                                  
#@description:  Generate relaxed candidate structures                                                       
#@author:     	Rebecca F. Alford                   
#@email:    		rfalford12@gmail.com                                          
###################################################################

## TODO: Rename sampled_candidates files

import random, os, sys
import hpc_util, read_config
from string import Template

def run_refine_decoys_calc( energy_fxn, restore, config, test_name, resolution, xml ): 
	"""
		A function for refining canddiate structures for the decoy discrimination test

		Arguments: 
			energy_fxn = energy function to use for calculations (typically, name of the weights file)
			config = path to benchmark, Rosetta executables
			targets = list of targets
			test_name = Name of test
	"""

	print( "Generating data for decoy discrimination test...")

	# Read list of test case IDs and PDBs
	targets_path = config.benchmark_path + "targets/structure/D6_decoy_discrimination/" + resolution + "/"
	list_of_targets = targets_path + "targets.list"

	with open( list_of_targets, 'rt' ) as f: 
			targets = f.readlines()
			targets = [ x.strip() for x in targets ]

	# Generate path to executable
	executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
	xml_script =  config.benchmark_path + "src/xml/" + xml

	# Change directories to a data analysis dir
	outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name + "/" + resolution + "/"  
	if ( not os.path.isdir( outdir ) ): 
		os.system( "mkdir " + outdir )
		os.chdir( outdir )

	# Iterate through each target
	for target in targets:

		print("Submitting refinement calculations for case:", target)

		# Read the list of decoy lists
		list_of_decoy_lists = []
		if ( resolution == "hires" ):
			list_of_decoy_lists = targets_path + "/" + target + "/list.of.decoy.lists"
		else: 
			list_of_decoy_lists = targets_path + "/" + target + "/" + target + "_model_subset_lists"

		with open( list_of_decoy_lists, 'rt' ) as f: 
			list_of_lists = f.readlines()
			list_of_lists = [ x.strip() for x in list_of_lists ]

		# Read general target variables
		spanfile = targets_path + "/" + target + "/" + target + ".span" 
		scorefile = target + "_refined.sc"
		casedir = outdir + "/" + target
		if ( not os.path.isdir( casedir ) ): 
			os.system( "mkdir " + casedir )
			os.chdir( casedir )

		# Iterate through each decoy list
		i = 0
		for dlist in list_of_lists: 

			i = i + 1

			# Setup case-specific variables
			models_list = targets_path + "/" + target + "/" + dlist

			# Generate a string of arguments from the case-specific variables
			s = Template( " -relax:constrain_relax_to_start_coords -in:file:l $models_list -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile $scorefile -out:path:all $outdir -nstruct $nmodels -run:multiple_processes_writing_to_one_directory ")
			arguments = s.substitute( models_list=models_list, span=spanfile, xml=xml_script, sfxn=energy_fxn, outdir=casedir, nmodels=5, scorefile=scorefile )

			# Restore/lipid composition flags
			if ( restore ): 
				arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior "
			else:
				arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:has_pore false "

			# Write jobfile and submit to the HPC
			jobname = target + "_refine_" + str(i)
			hpc_util.submit_condor_job( casedir, jobname, executable, arguments, 5 )


