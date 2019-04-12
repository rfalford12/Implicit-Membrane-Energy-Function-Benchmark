#!/usr/bin/env python
# @file: test_2.1_generate_flex_backbone_ddG.py
# @brief: Script for generating decoys for flexible backbone ddG calculation
# @notes: This script should **always** be run from the home membrane-efxn directory
# @author: Rebecca F. Alford (ralford3@jhu.edu)

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
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    os.system( "condor_submit " + filename )

def run_flexible_ddG_calc( energy_fxn, rosetta_exe_path, test_name, mutation_list, xml_protocol, restore, implicit_lipids ):
    """
    A general function for calculating the ddG of mutation with a flexible backbone and side chains

    Arguments: 
    - energy_fxn: energy function to use (typically, name of the weights file)
    - rosetta_exe_path = pah to compiled Rosetta executable
    - test_name = Name of energy landscpae test variant
    - mutations_list = List of input mutations
    - xml_protocol = Path to RosettaScript defining the landscape search protocol 
    - restore = Restore talaris behavior?
    - implicit_lipids = Use implicit lipid parameters?
    """
    
    print("Initializing flexible backbone ddG test for " + test_name)

    # Read list of mutations to generate
    inputs = benchmark + "inputs/test_2.1_ddG_of_mutation/" + test_name 
    list_of_mutations = inputs + "/" + mutation_list 
    with open( list_of_mutations, 'rb' ) as f: 
        mutations = f.readlines()
    mutations = [ x.strip() for x in mutations ]

    # Generate path to executable
    executable = rosetta_exe_path + "backrub" + "." + platform + compiler + buildenv

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/test_2.1_ddG_of_mutation"
    if ( not os.path.isdir( outdir ) ): 
        os.system( "mkdir " + outdir )
        os.chdir( outdir )

    outdir = outdir + "/" + test_name
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # Setup pdb and spanfiles
    pdbfile = inputs + "/" + test_name + "_min.pdb" 
    spanfile = inputs + "/" + test_name + ".span" 

    # For each test case, generate specific arguments, a condor file, and then run
    for mut in mutations:

        # Make a directory for the specific mutation
        mutationdir = outdir + "/" + mut
        os.mkdir( mutationdir )
        os.chdir( mutationdir )
        
        # Setup path to resfile
        resfile = inputs + "/resfiles/" + test_name + "_" + mut + ".resfile"
        scorefile = mutationdir + "/" + test_name + "_" + mut + "scores.sc"

        # Setup commandline arguments
        arguments = " -in:membrane -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -ex1 -ex2 -extrachi_cutoff 0 -backrub:ntrials 10000 -score:weights " + energy_fxn + " -nstruct 200 -run:multiple_processes_writing_to_one_directory -resfile " + resfile + " -out:file:scorefile " + scorefile + " -out:path:all " + mutationdir
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 20.0 -mp:lipids:composition DLPC"

        # Write arguments and executable to a separate file
        jobfile = mutationdir + "/" + mut + "_flex_bb_ddG.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate a condor submission file and submit the job to Jazz
        print("Submitting test case for " + test_name + ": " +  mut )
        write_and_submit_condor_script( mutationdir, mut, jobfile, "", 200 )

def main( args ):

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_t01 --which_tests all" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function weights file", )

    parser.add_option('--stable', '-s', 
        action="store", 
        help="For reference runs, use the stable Rosetta master branch", )

    parser.add_option( '--restore_talaris', '-r', 
        action="store", 
        help="Restore talaris behavior using tthe flag -restore_talaris_behavior for reference runs")

    parser.add_option( '--implicit_lipids', '-i', 
        action="store", 
        help="Use implicit lipids and default parameters when running this benchamrk", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options are provided and contents are valid
    if ( not Options.energy_fxn or not Options.stable ): 
        print("Missing required options --energy_fxn or --stable! Exiting...")
        sys.exit()

    # Choose rosetta path based on input options
    rosetta_exe_path = ""
    if ( Options.stable == "true" ): 
        rosetta_exe_path = rosettadir_stable
    else: 
        rosetta_exe_path = rosettadir

    # Set option for restoring talaris behaviorrun_fixed
    restore = False
    if ( Options.restore_talaris == "true" ): 
        restore = True

    include_lipids = False
    if ( Options.implicit_lipids == "true" ): 
        include_lipids = True

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

    ## Run all three fleixble ddG prediction tests
    # Mutations at the hydrophobic core in OmpLA (20 mutations, DLPC at RT)
    #run_flexible_ddG_calc( Options.energy_fxn, rosetta_exe_path, "OmpLA", "OmpLA_mutations.list", "xml/monomer_flex_ddG.xml", restore, include_lipids )
    
    # Mutations at the hydrophobic core in PagP (20 mutations, DLPC at RT)
    #run_flexible_ddG_calc( Options.energy_fxn, rosetta_exe_path, "PagP", "PagP_mutations.list", "xml/monomer_flex_ddG.xml", restore, include_lipids )
    
    # Mutations at varying depths to test aromatic preferences in OmpLA (38 mutations, DLPC at RT)
    #run_flexible_ddG_calc( Options.energy_fxn, rosetta_exe_path, "OmpLA_aro", "OmpLA_aro_mutations.list", "xml/monomer_flex_ddG.xml", restore, include_lipids )

    # Alanine reference mutations at varying positions
    run_flexible_ddG_calc( Options.energy_fxn, rosetta_exe_path, "OmpLA_aro_ref", "OmpLA_aro_ref_mutations.list", "xml/monomer_flex_ddG.xml", restore, include_lipids )

if __name__ == "__main__" : main(sys.argv)
