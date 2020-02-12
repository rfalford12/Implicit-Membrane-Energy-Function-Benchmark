#!/usr/bin/env python
###################################################################
#@file:         make_peptide_energy_landscape.py                                                                                   
#@description:  Generate orientation-dependent peptide energy map                                               
#@args:			--energy_fxn (wts) --targets (list of targets)                                                          
#@author: 		Rebecca F. Alford                   
#@email: 		rfalford12@gmail.com                                          
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_peptide_energy_landscape_calc( energy_fxn, config, test_name, targets, xml, pH_mode=False, pH=7 ): 

	"""
	A function to generate orientation-dependent peptide energy landscapes

	Arguments: 
		energy_fxn = energy function to use for calculations (typically, name of the weights file)
		config = path to benchmark, Rosetta executables
		test_name = name of test directory
		targets = Name of targets subdirectory (e.g. A4_designed_surface_ahelices)
		xml = Path to RosettaScript defining the protein energy landscape search protocol
	"""

	print( "Generating data for peptide tilt and rotation angle test..." ) 

	# Read list of energy landscape test cases
	targets_path = config.benchmark_path + "targets/" + targets + "/"
	list_of_targets = targets_path + "targets.list"
	with open( list_of_targets, 'rt' ) as f: 
		test_cases = f.readlines()
	test_cases = [ x.strip() for x in test_cases ]

	executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
	xml_script =  config.benchmark_path + "/" + xml

	# Change directories to a data analysis dir
	outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name
	if ( not os.path.isdir( outdir ) ): 
		os.system( "mkdir " + outdir )
		os.chdir( outdir )

	# For each test case, generate specific arguments, a condor file, and then run
	for case in test_cases:

		# Setup case-specific variables (pdbfile, spanfile, xmlargs)
		pdbfile = targets_path + "/" + case + "/" + case + ".pdb"
		spanfile = "single_TM_mode"

		# Default composition is DLPC
		arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script + " -mp:lipids:composition DLPC"
		if ( pH_mode ): 
			arguments = arguments + " -pH_mode true -value_pH " + str(pH)

		casedir = outdir + "/" + case
		if ( pH_mode ): 
			casedir = casedir + "_" + str(pH)
		if ( not os.path.isdir( casedir ) ): 
			os.system( "mkdir " + casedir )
			os.chdir( casedir )

		# Write jobfile and submit to the HPC
		print("Submitting orientation sampling calculations for energy landscape case:", case ) 
		jobname = case + "_protein_energy_landscape"
		jobfile = hpc_util.make_jobfile( casedir, case, executable, arguments )
		hpc_util.submit_condor_job( casedir, jobname, jobfile, "" )
