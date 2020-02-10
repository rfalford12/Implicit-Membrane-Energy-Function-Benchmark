#!/usr/bin/env python
###################################################################
#@file:         make_protein_energy_landscape.py                                                                                   
#@description:  Generate orientation-dependent peptide energy map                                                                                                       
#@author: 		Rebecca F. Alford                   
#@email: 		rfalford12@gmail.com                                          
###################################################################

import random
import sys, os
import hpc_util, read_config
from string import Template

def run_protein_energy_landscape_calc( energy_fxn, restore, config, test_name ): 
	"""
	A function to generate orientation-dependent peptide energy landscapes

	Arguments: 
		energy_fxn = energy function to use for calculations (typically, name of the weights file)
		config = path to benchmark, Rosetta executables
		test_name = name of test directory
		xml = Path to RosettaScript defining the protein energy landscape search protocol
	"""

	print( "Generating data for protein energy landscape test..." ) 

	# Read list of energy landscape test cases
	targets_path = config.benchmark_path + "targets/orientation/B1_multispan_proteins/"
	list_of_targets = targets_path + "targets.list"
	with open( list_of_targets, 'rt' ) as f: 
		test_cases = f.readlines()
	test_cases = [ x.strip() for x in test_cases ]

	# Generate path to executable
	executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
	xml_script =  config.benchmark_path + "/src/xml/make_depth_vs_tilt_energy_landscape.xml"

	# Change directories to a data analysis dir
	outdir = config.benchmark_path + "data/" + test_name + "/" + energy_fxn
	if ( not os.path.isdir( outdir ) ): 
		os.system( "mkdir " + outdir )
		os.chdir( outdir )

	# For each test case, generate specific arguments, a condor file, and then run
	for case in test_cases:

		# Setup case-specific variables (pdbfile, spanfile, xmlargs)
		pdbfile = targets_path + "/" + case + "/" + case + "_ignorechain.pdb"
		spanfile = "from_structure"

		# Default composition is DLPC
		arguments = " -overwrite -out:nooutput -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script 
		if ( restore ): 
			arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior "
		else:
			arguments = arguments + " -mp:lipids:composition DLPC"

		casedir = outdir + "/" + case
		if ( not os.path.isdir( casedir ) ): 
		  os.system( "mkdir " + casedir )
		  os.chdir( casedir )

		# Write jobfile and submit to the HPC
		print("Submitting protein orientation sampling calculations for case:", case ) 
		jobname = case + "_energy_landscape"
		jobfile = hpc_util.make_jobfile( casedir, case, executable, arguments )
		hpc_util.submit_condor_job( casedir, jobname, jobfile, "", 1 )