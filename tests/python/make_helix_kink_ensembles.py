#!/usr/bin/env python
###################################################################
#@file:         make_helix_kink_ensembles.py                                 
#@description:  Generate NMA ensembles of kinked structures                                                       
#@author:     	Rebecca F. Alford                   
#@email:    	rfalford12@gmail.com                                          
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template

def run_nma_ensemble_calc( config, xml ): 
	"""
		A function for generating normal model structural ensembles

		Arguments: 
			config = path to benchmark, Rosetta executables
			energy_fxn = energy function to use for calculations (typically, name of the weights file)
			xml = Name of XML file for running NMA
	"""

	print( "Generating data for helix kink test...")

	energy_fxn = "ref2015"
	test_name = "helix-kinks"

	# Read list of test case IDs and PDBs
	targets_path = config.benchmark_path + "targets/structure/D5_helix_kinks/"
	list_of_targets = targets_path + "targets.list"

	with open( list_of_targets, 'rt' ) as f: 
			targets = f.readlines()
			targets = [ x.strip() for x in targets ]

	# Generate path to executable
	executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
	xml_script =  config.benchmark_path + "src/xml/" + xml

	# Change directories to a data analysis dir
	outdir = config.benchmark_path + "data/" + test_name + "/" + energy_fxn + "/" 
	if ( not os.path.isdir( outdir ) ): 
		os.system( "mkdir " + outdir )
		os.chdir( outdir )

	# Iterate through each target
	for target in targets:

		# Read general target variables
		input_pdb = targets_path + target + "/" + target + ".pdb"
		scorefile = target + "_NMA_ensemble.sc"
		nmodels = 100
		casedir = outdir + "/" + target
		if ( not os.path.isdir( casedir ) ): 
			os.system( "mkdir " + casedir )
			os.chdir( casedir )

		# Generate a string of arguments from the case-specific variables
		s = Template( "-in:file:s $in_pdb -parser:protocol $xml -out:file:scorefile $scorefile -out:path:all $outdir -nstruct $nmodels -run:multiple_processes_writing_to_one_directory")
		arguments = s.substitute( in_pdb=input_pdb, xml=xml_script, sfxn=energy_fxn, outdir=casedir, nmodels=nmodels, scorefile=scorefile )

		# Add argument for local refine if applicable
		suffix = "_NMA_ensemble.sh"
		jobname = target + "_NMA_ensemble"

		# Write jobfile and submit to the HPC
		print("Submitting NMA calculations for case:", target)
		jobfile = hpc_util.make_jobfile( casedir, target, executable, arguments, suffix )
		hpc_util.submit_condor_job( casedir, jobname, jobfile, "", 100 )
