#!/usr/bin/env python 
# @file: run_benchmark_4_protein_design.py
# @brief: Master script for lipid composition dependent protein design
# @notes: This script should **always** be run from the membrane_efxn_benchamrk directory
# @author: Rebecca F. Alford (ralford3@jhu.edu)

## Example: python python/run_benchmark_4_protein_design.py --energy_fxn franklin2018 --stable false --restore false --lipid_composition DOPC --temperature 37 --subset hsapiens


import sys, os
import random
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

##############################################################################
### Global data for benchmark runs on Jazz
benchmark = "/home/ralford/membrane_efxn_benchmark/"
rosettadir = "/home/ralford/apps/Rosetta/main/source/bin/"
rosettadir_stable = "/home/ralford/apps/Rosetta-stable/main/source/bin/"
platform = "linux"
buildenv = "release"
compiler = "gcc"
##############################################################################

def write_and_submit_condor_script( path, name, executable, arguments, queue_no=1, high_mem=False ):

    # Given a test path and
    filename = path + "/" + name + ".condor"
    with open( filename, 'w' ) as f:
        f.write( "#!/bin/bash\n" )
        f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
        f.write( "Universe = vanilla\n" )
        f.write( "output = " + path + "/" + name + ".out\n" )
        f.write( "error = " + path + "/" + name + ".err\n" )
        if ( high_mem == True ): 
            f.write( "request_memory = 15000\n")
        f.write( "notify_user = rfalford12@gmail.com\n" )
        f.write( "Executable = " + executable + "\n" )
        f.write( "Arguments = " + arguments + "\n" )
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    os.system( "condor_submit " + filename )

def run_fixed_backbone_design_calc( energy_fxn, rosetta_exe_path, restore, lipid_composition, temperature, subset  ): 
    """
    A function for running the fixed backbone design calculations needed for the sequence recovery test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        restore = Restore to talalaris behavior and Lazaridis IMM1 behavior
		lipid_composition = name of the lipid composition to use
		temperature = Temperature corresponding to lipid compositon parameters
		subset = name of protein subset to evaluate
    """

    print( "Initializing fixed backbone design calculations for sequence recovery test" ) 

    # Read list of test cases
    inputs = benchmark + "inputs/test_3.1_sequence_recovery"
    list_of_test_cases = inputs + "/monomer_chains_" + subset + ".list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to the fixbb executable
    executable = rosetta_exe_path + "fixbb." + platform + compiler + buildenv

    # Make a directory for the subset and lipid composition
    outdir = benchmark + "data/" + energy_fxn + "/" + subset + "_" + lipid_composition + "_" + str(temperature)
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, condor files, and then run
    for case in test_cases:

        # Make one directory per case
        casedir = outdir + "/" + case
        os.system( "mkdir " + casedir )
        os.chdir( casedir )

        # Setup arguments by substitution
        pdbfile = inputs + "/" + case + "/" + case + "_tr_ignorechain.pdb"
        spanfile = inputs + "/" + case + "/" + case + "_tr.span"
        s = Template( " -in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res" )
        arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
        else: 
        	arguments = arguments + " -mp:lipids:composition " + lipid_composition + " -mp:lipids:temperature " + str(temperature) 
  
        # Write arguments and executable to a separate file
        jobfile = casedir + "/" + case + "_" + lipid_composition + "_" + str(temperature) + "_seqrecov.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate a condor submission file and submit the job to Jazz
        print("Submitting fixed backbone design calculation for sequence recovery case:", case) 
        queue_no = 1
        high_mem = True 
        write_and_submit_condor_script( casedir, case, executable, arguments, queue_no, high_mem )

def main( args ):

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_t01 --which_tests all" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function weights file", )

    parser.add_option('--stable', '-s', 
        action="store", 
        help="For reference runs, use the stable Rosetta master branch", )

    parser.add_option( '--restore', '-r', 
        action="store", 
        help="Restore talaris behavior using tthe flag -restore_talaris_behavior for reference runs")

    parser.add_option( '--lipid_composition', '-l', 
        action="store", 
        help="Lipid composition to use for design", )

    parser.add_option( '--temperature', '-a', 
        action="store", 
        help="Temperature for lipid composition calculation", )

    parser.add_option( '--subset', '-t', 
    	action="store", 
    	help="Subset of proteins to design, typically by species, taxonomy, or localization", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options are provided and contents are valid
    if ( not Options.energy_fxn or not Options.stable ): 
        print("Missing required options --energy_fxn and/or --cluster_type! Exiting...")
        sys.exit()

    # Choose rosetta path based on input options
    rosetta_exe_path = ""
    if ( Options.stable == "true" ): 
        rosetta_exe_path = rosettadir_stable
    else: 
        rosetta_exe_path = rosettadir

    # Set option for restoring talaris behaviorrun_fixed
    restore = False
    if ( Options.restore == "true" ): 
        restore = True

    # If it doesn't exist, make the output data directory
    datadir = benchmark + "data/"
    if ( not os.path.isdir(datadir) ): 
        os.system( "mkdir " + datadir )
    os.system( "cd " + datadir )

    # Make an output data directory for the specific energy function
    outdir = benchmark + "data/" + Options.energy_fxn
    if ( not os.path.isdir(outdir) ): 
        os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

    # Run the applicable design calculation
    run_fixed_backbone_design_calc( Options.energy_fxn, rosetta_exe_path, restore, Options.lipid_composition, Options.temperature, Options.subset )

if __name__ == "__main__" : main(sys.argv)

