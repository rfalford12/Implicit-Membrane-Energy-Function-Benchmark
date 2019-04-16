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
import make_asymm_docked_complexes, make_asymm_docked_complexes
import make_designed_protein_scaffolds, make_hybridized_kink_ensembles
import make_peptide_energy_landscape, make_protein_energy_landscape
import make_refined_decoys

from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def create_outdirs( energy_fxn, config ): 

	# If it doesn't exist, make the output data directory
    datadir = config.benchmark_path + "data/"
    if ( not os.path.isdir(datadir) ): 
        os.system( "mkdir " + datadir )
    os.system( "cd " + datadir )

    # Make an output data directory for the specific energy function
    outdir = datadir + energy_fxn
    if ( not os.path.isdir(outdir) ): 
        os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

def main( args ): 

	# Read options from the command line
	parser = OptionParser(usage="usage %prog --energy_fxn membrane_v0 --which_tests all --restore false" )
	parser.set_description(main.__doc__)

	parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function weights file", )

    parser.add_option( '--which_tests', '-w', 
        action="store", 
        help="Specify which test groups to run. Options are: ddG, landscape, prediction", )

    parser.add_option( '--restore_talaris', '-r', 
        action="store", 
        help="Restore talaris behavior using tthe flag -restore_talaris_behavior for reference runs")

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options have been provided
    if ( not Options.energy_fxn or not Options.which_tests ): 
    	print("Missing required options --energy_fxn and/or --which_tests" )
    	sys.exit()

	# Read path configuration file
	config = read_config.read_config()
	create_outdirs( config )

	# Check test categories
	all_categories = [ "orientation", "stability", "structure", "design" ]
	categories = []
	if ( Options.which_tests == "all" ): 
		categories = all_categories
	else: 
		types = Options.which_tests.split(",")
		for t in types: 
			categories.append( t )

	# Generate benchmark data for structure features tests
	if ( "structure" ):

		# Make native refined structures
		nstruct_refined = 20
		run_refinement_calc( Options.energy_fxn, config, "structure/D2_singe_pass_mp_complexes", "protein-protein-docking", "mp_relax.xml" )
		run_refinement_calc( Options.energy_fxn, config, "structure/D3_multi_pass_mp_complexes", "protein-protein-docking", "mp_relax.xml" )


        # Docking calculation for small homodimer set (Lomize et al. 2017)
        #run_refined_native_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "small-homodimer-set", restore, include_lipids, False )
        #run_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "small-homodimer-set", restore, include_lipids, False )
        #make_relaxed_homology_models( Options.energy_fxn, "small-homodimer-set", rosetta_exe_path, restore, include_lipids )


        # Docking calculation for large bound-bound set (Hurwitz et al. 2016)
        #run_refined_native_heterodimer_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-bound-set", restore, include_lipids, False )
        #run_heterodimer_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-bound-set", restore, include_lipids, False )
        #make_relaxed_homology_models( Options.energy_fxn, "large-bound-set", rosetta_exe_path, restore, include_lipids )

        # Docking calculation for large unbound set (simulated, Hurwitz et al. 2016)
        #make_relaxed_homology_models( Options.energy_fxn, "large-unbound-set", rosetta_exe_path, restore, include_lipids )
        #find_lowest_scoring_relaxed_native( Options.energy_fxn, "large-unbound-set" )
        #make_prepacked_docking_input( Options.energy_fxn, "large-unbound-set", rosetta_exe_path, restore, include_lipids )
        #run_heterodimer_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-unbound-set", restore, include_lipids, False )
    


if __name__ == "__main__" : main(sys.argv)