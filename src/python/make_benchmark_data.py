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

	# Read path configuration file
	config = read_config.read_config()
	create_outdirs( config )




if __name__ == "__main__" : main(sys.argv)