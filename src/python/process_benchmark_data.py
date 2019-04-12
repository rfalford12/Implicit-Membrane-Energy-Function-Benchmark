#!/usr/bin/env python
# @file: benchmark_analysis.py
# @brief: Master script for analyzing data generated by benchmark_analysis.py
# @notes: This script should **always** be run from the home membrane-efxn directory
# @author: Rebecca F. Alford (rfalford12@gmail.com)

import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

## Global data pertaining to the benchmark run ##
benchmark = "/home/ralford/membrane-efxn-benchmark/"
rosettadir = "/home/ralford/apps/Rosetta/main/source/bin/"
buildenv = "release"
compiler = "gcc"
boincdir = "/home/ralford/apps/Rosetta/bakerlab_scripts/boinc/"

def main( args ):

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_01--which_tests all" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function dataset to post process", )

    parser.add_option('--which_tests', '-w',
        action="store",
        help="Which tests should I analyze? all, intuition, or validation", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    if ( not Options.energy_fxn or not Options.which_tests ):
        print "Missing required options --energy_fxns and --which_tests! Exiting..."
        sys.exit()

    if ( Options.which_tests != 'all' and Options.which_tests != 'validation' ):
        print "Illegal value for option --which_tests! Found", Options.which_tests
        sys.exit()

    # Make an output data directory
    outdir = "/home/ralford/membrane-efxn-benchmark/analysis/"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    if ( Options.which_tests == 'all' or Options.which_tests == 'validation' ):
        print "Analyzing results from validation tests"
        #process_seq_recovery_results( Options.energy_fxn )
        #process_decoy_discrimination_results( Options.energy_fxn )
        process_docking_results( Options.energy_fxn )

def process_decoy_discrimination_results( energy_fxn ):

    print "Process decoy discrimination results generated with the energy function", energy_fxn

    # Read in list of alpha helical folding problems
    path_to_test = benchmark + "tests/validation/test-decoy-discrimination"
    list_of_test_cases = path_to_test + "/inputs/decoy_sets.list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Change directories to a data analysis dir
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/decoy-discrimination"
    os.chdir( outdir )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test_cases:

        print "Estimating decoy discrimination for", case

        outdir = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/decoy-discrimination/" + case
        os.chdir( outdir )

        # First, use Rosetta to rescore every single PDB relative to the native
        # Assuming a non-ref energy function for now
        # 

        # output_sc = outdir + "/" + case + "_rescored.sc"
        # os.system( "ls -d " + outdir + "/*_0001.pdb > " + models_list )
        # r = Template( " -overwrite -in:file:l $refined_list -in:file:native $native -in:membrane -mp:setup:spanfiles $spanfile -out:file:scorefile $output_sc -score:weights $energy_fxn" )
        # 
        # write_and_submit_condor_script( outdir, case + "_" + energy_fxn, executable, arguments0, 1 )

        # os.system( executable + arguments0 )

        # Specify paths to old and new scorefiles
        input_path = outdir + "/" + case + "_rescored.sc"
        output_path = outdir + "/" + case + "_rescored_refined.sc"

        # Parse each scorefile to exclude incomplete lines
        python_script = path_to_test + "/parse_scorefile.py"
        s = Template( " --num_elements 29 --input $input --output $output" )
        arguments1 = s.substitute( input=input_path, output=output_path )
        os.system( "python " + python_script + arguments1 )

        # Run decoy discrimination on each parsed scorefile (here, its by Irms and I_sc)
        scoring_script = boincdir + "score_energy_landscape.py"
        output_discfile = outdir + "/" + case + ".disc"
        t = Template( " -terms rms total_score -abinitio_scorefile $refined_sc > $refined_disc" )
        arguments2 = t.substitute( refined_sc=output_path, refined_disc=output_discfile )
        os.system( "python " + scoring_script + arguments2 )

def process_seq_recovery_results( energy_fxn ):

    print "Process sequence recovery results generated with the energy function", energy_fxn

    # Define a path to the test directory and native list
    path_to_test = benchmark + "tests/validation/test-seq-recovery"
    list_of_test_cases = path_to_test + "/inputs/monomer_chains_short.list"
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Make the native list from the existing list of monomer chains
    native_list = path_to_test + "/inputs/native_short.list"

    # Change directories into a data analysis dir
    outdir = "/home/ralford/membrane-efxn-benchmark/data/" + energy_fxn + "/seq-recovery"
    os.chdir( outdir )
    os.system( "ls -d $PWD/*/*_0001.pdb > redesigned.list" )
    redesign_list = outdir + "/redesigned.list"
    output_fn = outdir + "/" + energy_fxn + "_seqrecov.txt"

    # Setup sequence recovery analysis
    executable = rosettadir + "sequence_recovery.linux" + compiler + buildenv
    s = Template( " -overwrite -native_pdb_list $natives -redesign_pdb_list $redesigned -seq_recov_filename $output" )
    arguments = s.substitute( natives=native_list, redesigned=redesign_list, output=output_fn )
    write_and_submit_condor_script( outdir, "seqrecov", executable, arguments, 1 )

def write_and_submit_condor_script( path, name, executable, arguments, queue_no=1 ):

    # Given a test path and
    filename = path + "/" + name + ".condor"
    with open( filename, 'w' ) as f:
        f.write( "#!/bin/bash\n" )
        f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
        f.write( "Universe = vanilla\n" )
        f.write( "output = " + path + "/" + name + ".out\n" )
        f.write( "error = " + path + "/" + name + ".err\n" )
        f.write( "notify_user = rfalford12@gmail.com\n" )
        f.write( "Executable = " + executable + "\n" )
        f.write( "Arguments = " + arguments + "\n" )
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    os.system( "condor_submit " + filename )

if __name__ == "__main__" : main(sys.argv)
