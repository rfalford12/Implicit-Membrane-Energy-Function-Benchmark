#!/usr/bin/env python
""" Master script for generating benchmark data

This module generates data for any of the twelve implicit 
membrane energy function benchmark tests. This is the first
step of three for executing and evaluating the scientific tests.

Authors: 
	Rebecca Alford <ralford3@jhu.edu> 

Example: 
	$ python generate_benchmark_data.py --energy_fxn franklin2019
	--which_tests all --restore_talaris False

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- which_tests: Run all tests, or specify by comma-separated list

Requirements: 
	- Rosetta release 246 or greater
	- PyRosetta4 for Python 3.6 or 3.7
"""

import sys, os, random
import hpc_util, read_config
import make_helix_kink_ensembles
import make_peptide_energy_landscape
import make_protein_energy_landscape
import make_asymm_docked_complexes
import make_designed_protein_scaffolds
import make_refined_decoys
import predict_hydrophobic_length
import predict_ddG
import predict_side_chain_distribution

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

all_tests = [ "sc-distribution", "ddG-of-insertion", "ddG-of-mutation", "ddG-of-pH-insertion", "decoy-discrimination", "tm-peptide-tilt", "helix-kinks", "hydrophobic-length", "adsorbed-peptide-tilt-angle", "protein-protein-docking", "protein-tilt-angle", "sequence-recovery" ]

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
	if ( Options.energy_fxn == "franklin2019" or Options.energy_fxn == "ref2015" ): 
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

	# Test #1: Tilt angles for transmembrane peptides
	if ( "tm-peptide-tilt-angle" in test_names ): 

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "tm-peptide-tilt-angle", "tilt_angle/A1_native_tm_ahelices", "src/xml/make_depth_vs_tilt_energy_landscape.xml" )
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "tm-peptide-tilt-angle", "tilt_angle/A3_designed_tm_ahelices", "src/xml/make_depth_vs_tilt_energy_landscape.xml" )

	# Test #2: Rotation angles for surface adsorbed peptides
	if ( "adsorbed-peptide-tilt-angle" in test_names ): 
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "adsorbed-peptide-tilt-angle", "tilt_angle/A2_native_surface_ahelices", "src/xml/make_depth_vs_helix_rot_energy_landscape.xml" )

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "adsorbed-peptide-tilt-angle", "tilt_angle/A4_designed_surface_ahelices", "src/xml/make_depth_vs_helix_rot_energy_landscape.xml" )

	# Test #3: Orientation of membrane proteins with complex topologies
	if ( "protein-tilt-angle" in test_names ): 

		make_protein_energy_landscape.run_protein_energy_landscape_calc( Options.energy_fxn, restore, config, "protein-tilt-angle" )

	# Test #4: Membrane protein hydrophobic thickness
	if ( "hydrophobic-length" in test_names ): 

		predict_hydrophobic_length.run_hydrophobic_length_calc( Options.energy_fxn, restore, config, "hydrophobic-length" )

	# Test #5: Stability of transmembrane peptides at neutral pH
	if ( "ddG-of-insertion" in test_names ): 

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "ddG-of-pH-insertion", "stability/C4_polyLeu_helical_peptides", "src/xml/make_depth_vs_tilt_energy_landscape.xml" )

	# Test #6: Stability of transmembrane peptides at acidic pH
	if ( "ddG-of-pH-insertion" in test_names ): 

		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "ddG-of-pH-insertion", "stability/C5_pHLIP_helical_peptides", "src/xml/make_depth_vs_tilt_pH_energy_landscape.xml", True, 4 )
		make_peptide_energy_landscape.run_peptide_energy_landscape_calc( Options.energy_fxn, config, "ddG-of-pH-insertion", "stability/C5_pHLIP_helical_peptides", "src/xml/make_depth_vs_tilt_pH_energy_landscape.xml", True, 8 )

	# Test #7: ddG of mutation
	if ( "ddG-of-mutation" in test_names ): 

		predict_ddG.run_ddG_of_mutation_calc( config, Options.energy_fxn, "C1_OmpLA_canonical_ddGs", "1qd6.pdb", "1qd6.span", "OmpLA_Moon_Fleming_set.dat" )
		predict_ddG.run_ddG_of_mutation_calc( config, Options.energy_fxn, "C2_PagP_canonical_ddGs", "3gp6_A.pdb", "3gp6_A.span", "PagP_Marx_Fleming_set.dat" )
		predict_ddG.run_ddG_of_mutation_calc( config, Options.energy_fxn, "C3_OmpLA_aro_ddGs", "1qd6.pdb", "1qd6.span", "OmpLA_aro_McDonald_Fleming_set.dat" )		

	# Test #8 and #9: Sequence Recovery and side chain distribution
	if ( ("sequence-recovery" or "sc-distribution") in test_names ): 

		make_designed_protein_scaffolds.run_fixed_backbone_design_calc( config, Options.energy_fxn,  "targets", "sequence-recovery" )

	# Test #10: Native structure discrimination
	if ( "decoy-discrimination" in test_names ): 

		make_refined_decoys.run_refine_decoys_calc( Options.energy_fxn, restore, config, "decoy-discrimination", "hires", "mp_relax.xml" ) 
		make_refined_decoys.run_refine_decoys_calc( Options.energy_fxn, restore, config, "decoy-discrimination", "lowres", "mp_relax.xml" )

	# Test #11: Membrane helix kinks
	if ( "helix-kinks" in test_names ): 

		make_helix_kink_ensembles.run_nma_ensemble_calc( config, "nma.xml")

	# # Test #12: Protein-protein docking
	# if ( "protein-protein-docking" in test_names ): 

	# 	# Step 12.1 - Generate prepacked structures
	# 	make_asymm_docked_complexes.run_prepack_calc( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking" )
	# 	make_asymm_docked_complexes.run_prepack_calc( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking" )
	# 	make_asymm_docked_complexes.run_prepack_calc( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking" )

	# 	# Step 12.2 - Remove "MEM" from the PDB
	# 	make_asymm_docked_complexes.post_process_prepack_pdb( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True )
	# 	make_asymm_docked_complexes.post_process_prepack_pdb( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking" )
	# 	make_asymm_docked_complexes.post_process_prepack_pdb( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking" )

	# 	# Step 12.3 - Generate 5K docking decoys per target
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", False, True )
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", False  )
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", False  )

	# 	# Step 12.4 - Generate 100 refined native docking decoys per target
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True, True )
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", True  )
	# 	make_asymm_docked_complexes.run_docking_calc( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", True  )



if __name__ == "__main__" : main(sys.argv)
