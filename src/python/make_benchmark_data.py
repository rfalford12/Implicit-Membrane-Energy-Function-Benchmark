#!/usr/bin/env python
###################################################################
#@file:         make_benchmark_data.py                                                                                  
#@description:  Generate data required for membrane efxn benchmark suite                                              
#@args:			--energy_fxn (wts)                                                
#@author: 		Rebecca F. Alford                   
#@email: 		rfalford12@gmail.com                                          
###################################################################

import sys, os, random
import hpc_util, read_config
import make_refined_decoys

## TODO for later
#import make_asymm_docked_complexes
#import make_designed_protein_scaffolds, make_hybridized_kink_ensembles
#import make_peptide_energy_landscape, make_protein_energy_landscape

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

all_tests = [ "aro-distribution", "charge-distribution", "ddG-of-insertion", "ddG-of-mutation", "ddG-of-pH-insertion", "decoy-discrimination", "energy-function-landscape", "helix-kinks", "hydrophobic-length", "peptide-tilt-angle", "protein-protein-dokcing", "protein-tilt-angle", "sequence-recovery", "symmetric-protein-docking" ]

def create_outdirs( energy_fxn, config, list_of_tests ): 

	# If it doesn't exist, make the output data directory
	datadir = config.benchmark_path + "data/"
	if ( not os.path.isdir(datadir) ): 
		os.system( "mkdir " + datadir )
	os.system( "cd " + datadir )

	# Create test data directories
	for test in list_of_tests: 
		testdir = datadir + test + "/"
		print(testdir)
		if ( not os.path.isdir(testdir) ): 
			os.system( "mkdir " + testdir )

		efxndir = testdir + energy_fxn
		print(efxndir)
		if ( not os.path.isdir(efxndir) ): 
			os.system( "mkdir " + efxndir )

def main( args ): 

	# Read options from the command line
	parser = OptionParser(usage="usage %prog --energy_fxn membrane_v0 --which_tests all --restore false" )
	parser.set_description(main.__doc__)
	
	parser.add_option( '--energy_fxn', '-e', action="store", help="Name of energy function weights file", )
	parser.add_option( '--which_tests', '-w', action="store", help="Specify tests run (comma separated list)", )
	parser.add_option( '--restore_talaris', '-r', action="store", help="Restore talaris behavior using tthe flag -restore_talaris_behavior for reference runs", )

	(options, args) = parser.parse_args(args=args[1:])
	global Options
	Options = options

	# Check that required options have been provided
	if ( not Options.energy_fxn or not Options.which_tests ): 
		print("Missing required options --energy_fxn and/or --which_tests" )
		sys.exit()

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
				sys.exit( "No such test " + name + ". Available tests are " + all_tests + ". Exiting!" )

	# Create directories based on test names and user-specified energy function
	create_outdirs( Options.energy_fxn, config, test_names )

	# Generate benchmark data for structure features tests
	# if ( "structure" ):

	# 	##### Test: Protein-protein docking ########

	# 	# Make native refined structures
	# 	nstruct_native_refined = 20
	# 	make_refined_decoys.run_refinement_calc( Options.energy_fxn, config, "structure/D2_singe_pass_mp_complexes", "protein-protein-docking/", "mp_relax.xml", nstruct_native_refined, "D2_single_tm" )
	# 	make_refined_decoys.run_refinement_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes/bound-complexes", "protein-protein-docking", "mp_relax.xml", nstruct_native_refined, "D3_bound" )
	# 	make_refined_decoys.run_refinement_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes/unbound-complexes", "protein-protein-docking", "mp_relax.xml", nstruct_native_refined, "D3_unbound" )

	# 	# Make protein-protein complex models via rigid body docking
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "structure/D2_singe_pass_mp_complexes", "protein-protein-docking", "D2_single_tm" )
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes/bound-complexes/", "protein-protein-docking", "bound" )
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes/unbound-complexes/", "protein-protein-docking", "unbound" )


if __name__ == "__main__" : main(sys.argv)
