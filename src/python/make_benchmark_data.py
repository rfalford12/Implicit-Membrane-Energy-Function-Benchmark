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
import make_helix_kink_ensembles
import make_peptide_energy_landscape
import make_protein_energy_landscape
import make_asymm_docked_complexes
import make_refined_decoys
import predict_hydrophobic_length
import predict_ddG
import predict_side_chain_distributions

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

all_tests = [ "sc-distribution", "ddG-of-insertion", "ddG-of-mutation", "ddG-of-pH-insertion", "decoy-discrimination", "energy-function-landscape", "helix-kinks", "hydrophobic-length", "peptide-tilt-angle", "protein-protein-docking", "protein-tilt-angle", "sequence-recovery", "symmetric-protein-docking" ]

def create_outdirs( energy_fxn, config, list_of_tests ): 

	# If it doesn't exist, make the output data directory
	datadir = config.benchmark_path + "data/"
	if ( not os.path.isdir(datadir) ): 
		os.system( "mkdir " + datadir )
	os.system( "cd " + datadir )

	# Create test data directories
	for test in list_of_tests: 
		testdir = datadir + test + "/"
		if ( not os.path.isdir(testdir) ): 
			os.system( "mkdir " + testdir )

		efxndir = testdir + energy_fxn
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

	# Set restore variable based on energy function type
	restore = True
	if ( Options.energy_fxn != "franklin2019" or Options.energy_fxn != "ref2015" ): 
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

	# Create directories based on test names and user-specified energy function
	create_outdirs( Options.energy_fxn, config, test_names )

	# Test ##: Generate benchamrk data for hydrophobic length test
	if ( "hydrophobic-length" in test_names ): 

		predict_hydrophobic_length.run_hydrophobic_length_calc( Options.energy_fxn, restore, config, "hydrophobic-length" )

	# Test ##: Generate benchamrk data for peptide tilt and rotation angle test
	if ( "peptide-tilt-angle" in test_names ): 

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "peptide-tilt-angle", "tilt_angle/A1_native_tm_ahelices", "src/xml/make_depth_vs_tilt_energy_landscape.xml" )
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "peptide-tilt-angle", "tilt_angle/A2_native_surface_ahelices", "src/xml/make_depth_vs_helix_rot_energy_landscape.xml" )
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "peptide-tilt-angle", "tilt_angle/A3_designed_tm_ahelices", "src/xml/make_depth_vs_tilt_energy_landscape.xml" )
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "peptide-tilt-angle", "tilt_angle/A4_designed_surface_ahelices", "src/xml/make_depth_vs_helix_rot_energy_landscape.xml" )

	# Test ##: Generate benchamrk data for protein energy landscape test
	if ( "protein-tilt-angle" in test_names ): 

		make_protein_energy_landscape.run_protein_energy_landscape_calc( Options.energy_fxn, restore, config, "protein-tilt-angle" )

	# Test ##: Generate benchamrk data for pH-sensitive insertion
	if ( "ddG-of-pH-insertion" in test_names ): 

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "ddG-of-pH-insertion", "stability/C5_pHLIP_helical_peptides", "src/xml/make_depth_vs_tilt_pH_energy_landscape.xml", True, 4 )
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "ddG-of-pH-insertion", "stability/C5_pHLIP_helical_peptides", "src/xml/make_depth_vs_tilt_pH_energy_landscape.xml", True, 8 )

	# Test ##: Generate benchamrk data for constant pH insertion
	if ( "ddG-of-insertion" in test_names ): 

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "ddG-of-pH-insertion", "stability/C4_polyLeu_helical_peptides", "src/xml/make_depth_vs_tilt_energy_landscape.xml" )

	# Test ##: Generate benchamrk data for ddG of mutation predictions
	if ( "ddG-of-mutation" in test_names ): 

		predict_ddG.run_ddG_of_mutation_calc( config, Options.energy_fxn, "C1_OmpLA_canonical_ddGs", "1qd6.pdb", "1qd6.span", "OmpLA_Moon_Fleming_set.dat" )
		predict_ddG.run_ddG_of_mutation_calc( config, Options.energy_fxn, "C2_PagP_canonical_ddGs", "3gp6_A.pdb", "3gp6_A.span", "PagP_Marx_Fleming_set.dat" )
		predict_ddG.run_ddG_of_mutation_calc( config, Options.energy_fxn, "C3_OmpLA_aro_ddGs", "1qd6.pdb", "1qd6.span", "OmpLA_aro_McDonald_Fleming_set.dat" )		

	# Test ##: Generate benchmark data for decoy-discrimination tests
	if ( "decoy-discrimination" in test_names ): 

		make_refined_decoys.run_refine_decoys_calc( Options.energy_fxn, restore, config, "decoy-discrimination", "hires", "mp_relax.xml" ) 
		make_refined_decoys.run_refine_decoys_calc( Options.energy_fxn, restore, config, "decoy-discrimination", "lowres", "mp_relax.xml" )

	# Test ##: Generate benchmark data for docking test
	if ( "helix-kinks" in test_names ): 

		make_helix_kink_ensembles.run_nma_ensemble_calc( config, "nma.xml")
        
        # Test ##: Generate benchmark data for side chain distribution test
        if ( "sc-distribution" in test_names ): 
                path = "/home/ralford/Implicit-Membrane-Energy-Function-Benchmark/data/sequence-recovery/franklin2019"
                predict_side_chain_distribution.compute_side_chain_distribution( config, path + "natives.list", path + "redesigned_allpath.list" )

	# Test ##: Generate benchmark data for docking test
	if ( "protein-protein-docking" in test_names ): 

		# Substep ##1 - Generate prepacked structures
		#make_asymm_docked_complexes.run_prepack_calc( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking" )
		#make_asymm_docked_complexes.run_prepack_calc( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking" )
		#make_asymm_docked_complexes.run_prepack_calc( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking" )

		# Substep ##2 - Remove "MEM" from the PDB
		#make_asymm_docked_complexes.post_process_prepack_pdb( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True )
		#make_asymm_docked_complexes.post_process_prepack_pdb( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking" )
		#make_asymm_docked_complexes.post_process_prepack_pdb( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking" )

		# Substep ##3 - Generate 5K docking decoys per target
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", False, True )
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", False  )
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", False  )

		# Substep ##4 - Generate 100 refined native docking decoys per target
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True, True )
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", True  )
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", True  )

<<<<<<< HEAD
=======
		# Make protein-protein complex models via rigid body docking
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "structure/D2_singe_pass_mp_complexes", "protein-protein-docking", "D2_single_tm" )
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes/bound-complexes/", "protein-protein-docking", "D3_bound" )
		make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes/unbound-complexes/", "protein-protein-docking", "D3_unbound" )
>>>>>>> 7a37ab5789a2dc125cba0c57befd1068ceb2aeba


if __name__ == "__main__" : main(sys.argv)
