#!/usr/bin/env python
""" Master script for processing sequence and structure benchmarks

This module analyzes data from tests that probe sequence and 
structure features. This is the second step of three for executing
and evaluating the scientific testds. 

Authors: 
	Rebecca Alford <ralford3@jhu.edu> 

Example: 
	$ python analyze_benchmark_data.py --energy_fxn franklin2019
	--which_tests all --restore_talaris False

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- which_tests: Run all tests, or specify by comma-separated list

Requirements: 
	- Rosetta release 246 or greater
	- PyRosetta4 for Python 3.6 or 3.7
	- KinkFinder
"""

import sys, os, random, read_config
import make_asymm_docked_complexes
import make_designed_protein_scaffolds
import make_refined_decoys
import predict_side_chain_distribution

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

all_tests = [ "sc-distribution", "ddG-of-insertion", "ddG-of-mutation", "ddG-of-pH-insertion", "decoy-discrimination", "tm-peptide-tiolt", "helix-kinks", "hydrophobic-length", "adsorbed-peptide-tilt-angle", "protein-protein-docking", "protein-tilt-angle", "sequence-recovery" ]
sequence_tests = [ "sc-distribution", "sequence-recovery" ]
structure_tests = [ "decoy-discrimination", "helix-kinks", "protein-protein-docking" ]

def main( args ): 

	# Read options from the command line
	parser = OptionParser(usage="usage %prog --energy_fxn franklin2019 --which_tests all --restore_talaris false" )
	parser.set_description(main.__doc__)
	
	parser.add_option( '--energy_fxn', '-e', action="store", help="Name of energy function weights file", )
	parser.add_option( '--which_tests', '-w', action="store", help="Specify tests run (comma separated list)", )

	(options, args) = parser.parse_args(args=args[1:])
	global Options
	Options = options

	# Check that required options have been provided
	if ( not Options.energy_fxn or not Options.which_tests ): 
		print("Missing required options --energy_fxn and/or --which_tests" )
		sys.exit()

	# Set restore variable based on energy function type
	restore = True
	if ( Options.energy_fxn == "franklin2019" or Options.energy_fxn == "ref2015" or Options.energy_fxn "ref2015_memb" ): 
		restore = False 

	# Read path configuration file
	config = read_config.read_config()

	# Check test categories
	test_names = []
	if ( Options.which_tests == "all" ): 
		test_names = all_tests
	else: 
		test_names = Options.which_tests.split(",")
		# check that all names are valid
		for name in test_names: 
			if name not in all_tests: 
				sys.exit( "No such test " + name + ". Exiting!" )

	# Begin analysis steps for benchmarks 
	# this is a giant TODO that needs to be automated for many purposes :)

	# Test #9: Side chain distribution calculations
	# this is going to become an analysis step... 
	if ( "sc-distribution" in test_names ): 

		path = "/home/ralford/Implicit-Membrane-Energy-Function-Benchmark/data/sequence-recovery/franklin2019"
		predict_side_chain_distribution.compute_side_chain_distribution( config, path + "/natives.list", path + "/redesigned_allpath.list" )

	# protein-protein docking interfae analysis steps

			# Step 12.5 - Analyze docking decoys
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", False, True )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", False  )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", False  )

		# Step 12.6 - Analyze local refine decoys
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True, True )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", True  )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", True  )

	# for decoy discrimination, automate score_energy_landscape step and extrapolate results

	# sequence recovery - measure sequence recovery step

