#!/usr/bin/env python
""" Predict orientation of helical peptides

This script runs the MembraneEnergyLandscapeSampler on each target
in the dataset. The result is a complete mapping of energies to 
peptide orientation, given as a function of depth and tilt 
either at rotation angle=0 or minimized over all rotation 
angles. This script generates data for Tests 1,2,5,6.

Authors: 
	Rebecca Alford <ralford3@jhu.edu> 
	Rituparna Samanta <rituparna@utexas.edu>
Example: 
	$ import make_peptide_energy_landscape
	$ make_peptide_energy_landscape.run_peptide_energy_landscape_calc( 
	  energy_fxn, config, test_name, targets_dir, xml )

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- config: Container with path to benchmark and rosetta files
	- test_name: Name of benchmark test
	- targets_dir: Location of test cases
	- xml: Path to MembraneEnergyLandscapeSampler XML application

Requirements: 
	- Rosetta release 246 or greater
"""

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

	executable = config.rosetta_path + "/rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
	xml_script =  config.benchmark_path + xml
	print("xml_scrip:~",xml_script)

	# Change directories to a data analysis dir
	outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name
	print("outdir:~",outdir)

	if ( not os.path.isdir( outdir ) ): 
		os.system( "mkdir " + outdir )
		os.chdir( outdir )

	# For each test case, generate specific arguments, a condor file, and then run
	for case in test_cases:

		# Setup case-specific variables (pdbfile, spanfile, xmlargs)
		pdbfile = targets_path + case + "/" + case + ".pdb"
		spanfile = "single_TM_mode"
		print("pdbfile:~",pdbfile)
		
		#flag_axis=2:eigenvector; 0:vector joining top and bottom centers; 1:average of all tm-axis
		flag_axis = ["0.0"]

		# Default composition is DLPC
		if ( "tm-peptide-tilt-angle" in test_name ):
			interface = "0"
		if ( "adsorbed-peptide-tilt-angle" in test_name ):
			interface = "1"
			if( "tilt_angle/A2_native_surface_ahelices" in targets ):
				flag_axis = ["2.0"]

		if ( "ddG-of-insertion" in test_name ):
			interface = "0"
		if ( "ddG-of-pH-insertion" in test_name ):
			interface = "0"
		
		
		start_z = ["-60", "-40", "-20", "0", "20", "40"]
		end_z = ["-40", "-20", "0", "20", "40", "60"]

		for j in range(len(start_z)):
 
			#the default lipid composition is DLPC, based on the test, the lipid composition and temperature needs to be changed. 
			arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:script_vars interface=" + interface + " -parser::script_vars start_z=" + start_z[j] + " -parser::script_vars end_z=" + end_z[j] + " -parser::script_vars flag_axis=" + flag_axis[0] + " -parser:protocol " + xml_script + " -mp:lipids:composition DLPC -mp:lipids:temperature 20.0"
		
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
			jobname = case + "_" + end_z[j] + "_" + start_z[j] + "_protein_energy_landscape"
			print( jobname) 
			jobfile = hpc_util.make_jobfile( casedir, jobname, executable, arguments )
			#print( "outdir:~", outdir )
			print( "casedir:~", casedir )
			print( "case:~", case )
			#print( "executable:~", executable )
			#print( "arguments:~", arguments )
			#hpc_util.submit_condor_job( casedir, jobname, jobfile, "" )
