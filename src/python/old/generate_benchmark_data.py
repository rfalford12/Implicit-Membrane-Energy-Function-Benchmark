#!/usr/bin/env python
# @file: generate_benchmark_data_marcc.py
# @brief: Master script for running membrane energy function benchmarks on MARCC
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
### Global data for benchmark runs on MARCC
# benchmark = "/home-4/ralford3@jhu.edu/work/ralford3@jhu.edu/membrane_efxn_benchmark/"
# rosettadir = "/home-4/ralford3@jhu.edu/work/ralford3@jhu.edu/Rosetta/main/source/bin/"
# rosettadir_stable = "/home-4/ralford3@jhu.edu/work/ralford3@jhu.edu/Rosetta-stable/main/source/bin/"
# platform = "mpi.linux" 
# buildenv = "release"
# compiler = "gcc"
##############################################################################
### Global data for benchmark runs on XSEDE-stampede2
# benchmark = "/work/04819/ralford3/membrane_efxn_benchmark/"
# rosettadir = "/work/04819/ralford3/Rosetta-stable/main/source/bin/"
# rosettadir_stable = "/work/04819/ralford3/Rosetta-stable/main/source/bin/"
# platform = "mpi.linux" 
# buildenv = "release"
# compiler = "icc"
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

def write_and_submit_stampede2_slurm_script( path, name, jobfile, num_nodes=1 ): 
  
    # Create a new sbatch file named for the job type and test
    filename = path + "/" + name + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Stampede2 with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH -J " + name + "\n" )
        f.write( "#SBATCH -p normal\n" )
        f.write( "#SBATCH -N " + str(num_nodes) + "\n" )
        f.write( "#SBATCH -n 16\n" )
        f.write( "#SBATCH -t 24:0:0\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH -o " + path + "/" + name + ".%j.out\n" )
        f.write( "#SBATCH -e " + path + "/" + name + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rfalford12@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "#SBATCH -A TG-MCB180056\n" )
        f.write( "\n" )

        # Specify required modules
        f.write( "module load intel/18.0.0\n" )
        f.write( "export MKL_MIC_ENABLE=1\n" )

        # Provide a description of the job
        f.write(  "echo Starting MPI job running " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpiexec bash " + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    sbatch_command = "sbatch " + filename
    os.system( sbatch_command )


def write_and_submit_slurm_batch_script( path, name, jobfile, num_nodes=1 ): 

    # Create a new sbatch file named for the job type and test
    filename = path + "/" + name + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on MARCC with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + name + "\n" )
        f.write( "#SBATCH --partition=parallel\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --time=60:0:0\n" )
        f.write( "#SBATCH --mem=120GB\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + name + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + name + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rfalford12@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        f.write( "module unload openmpi gcc\n" )
        f.write( "module load intel-mpi git\n" )
        f.write( "ml --gcc\n")
        f.write( "ml gcc/4.8.2\n")

        # Provide a description of the job
        f.write(  "echo Starting MPI job running " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpirun ./" + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    sbatch_command = "sbatch " + filename
    os.system( sbatch_command )

def run_energy_landscape_calc( energy_fxn, rosetta_exe_path, cluster_type, test_name, input_list, xml_protocol, restore, implicit_lipids, single_TM="false", pH="0" ): 
    """
    A general functions for running energy landscape calculations given a list of input helices

    Arguments:
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        test_name = Name of energy landscpae test variant
        input_list = List of input helices
        xml_protocol = Path to RosettaScript defining the landscape search protocol
        restore = Restore talaris behavior? 
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
        single_TM = use the single_TM_peptide spanfile option (default false)
        pH = run the simulation at the user specified pH value (defualt 0)
    """

    print("Initializing energy landscape test for " + test_name)

    # Read list of energy landscape test cases
    inputs = benchmark + "inputs/" + test_name 
    list_of_test_cases = inputs + "/" + input_list 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate path to executable
    executable = rosetta_exe_path + "rosetta_scripts" + "." + platform + compiler + buildenv
    xml_script =  benchmark + xml_protocol

    # Change directories to a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/" + test_name
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = inputs + "/" + case + "/" + case + ".pdb"

        # Should I use a single_TM_peptide span estimation or a user-provided spanfile? 
        spanfile = ""
        if ( single_TM == "false" ): 
            spanfile = inputs + "/" + case + "/" + case + ".span"
        else: 
            spanfile = "single_TM_mode"

        # Should I tune the pH of my simulation?
        arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script
        if ( pH != "0" ): 
            arguments = arguments + " -pH_mode true -value_pH " + pH
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"

        # Write arguments and executable to a separate file
        jobfile = outdir + "/" + case + "_energy_landscape.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate a condor submission file and submit the job to Jazz
        print("Submitting test case for " + test_name + ": " +  case)
        if ( cluster_type == "MARCC" ): 
            write_and_submit_slurm_batch_script( outdir, case, jobfile )
        elif ( cluster_type == "STAMPEDE" ): 
            write_and_submit_stampede2_slurm_script( outdir, case, jobfile )
        else: 
            write_and_submit_condor_script( outdir, case, executable, arguments )

def run_ddG_of_mutation_calc( energy_fxn, list_of_ddGs, test_name, restore, implicit_lipids, aqueous_pore ): 
    """
    A general function for calculating the ddG of single point mutations

    Arguments:  
        energy_fxn = energy function to run the protocol with
        list_of_ddGs = Relative path to list of ddGs to calculate
        test_name = Name of benchmark set for calculating ddGs of mutation
    """

    print("Initializing ddG-of-mutation test for set " + test_name)

    # Specify path to test, setup new directories
    path_to_test = benchmark + "inputs/test_2.1_ddG_of_mutation"
    outdir = benchmark + "data/" + energy_fxn + "/test_2.1_ddG_of_mutation"
    os.system( "mkdir " + outdir )
    os.system( "cd " + outdir )

    # Specify path to python script and arguments
    python_script = benchmark + "/python/test_2.1_predict_ddG.py"
    energy_function = energy_fxn
    mlist = path_to_test + "/" + list_of_ddGs
    s = Template( "--energy_fxn $energy_func --mutation_list $list_of_mutations --outdir  $outdir --implicit_lipids $ilm ")
    arguments = s.substitute( energy_func=energy_function, list_of_mutations=mlist, outdir=outdir, ilm=implicit_lipids )
    if ( restore == "true" ): 
        arguments = arguments + " --restore True"

    print("Submitting ddG of mutation test case " + test_name)
    os.system( "python " + python_script + " " + arguments )

def run_fixed_backbone_design_calc( energy_fxn, rosetta_exe_path, cluster_type, restore, implicit_lipids, aqueous_pore  ): 
    """
    A function for running the fixed backbone design calculations needed for the sequence recovery test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        restore = Restore to talalaris behavior prior to ref2015 for reference benchmark run
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print("Initializing fixed backbone design calculations for sequence recovery test") 

    # Read list of test cases - typically full length transmembrane proteins
    inputs = benchmark + "inputs/test_3.1_sequence_recovery"
    list_of_test_cases = inputs + "/monomer_chains_ecoli.list" 
    with open( list_of_test_cases, 'rb' ) as f:
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate a path to the fixbb executable
    executable = rosetta_exe_path + "fixbb." + platform + compiler + buildenv

    # Change directories into a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/test_3.1_sequence_recovery"
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

    # For each test case, generate specific arguments, condor files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.1_sequence_recovery/" + case
        os.system( "mkdir " + outdir )
        os.chdir( outdir )

        # Setup arguments by substitution
        pdbfile = inputs + "/" + case + "/" + case + "_tr_ignorechain.pdb"
        spanfile = inputs + "/" + case + "/" + case + "_tr.span"
        s = Template( " -in:file:s $pdbfile -mp:setup:spanfiles $spanfile -score:weights $sfxn -in:membrane -out:path:all $outdir -in:file:load_PDB_components false -in:ignore_unrecognized_res" )
        arguments = s.substitute( pdbfile=pdbfile, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPE"
        if ( aqueous_pore == True ): 
            arguemnts = arguments + " -mp:pore:accomodate_pore true"

        # Write arguments and executable to a separate file
        jobfile = outdir + "/" + case + "_seqrecov.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate a condor submission file and submit the job to Jazz
        print("Submitting fixed backbone design calculation for sequence recovery case:", case)
        if ( cluster_type == "MARCC" ): 
            write_and_submit_slurm_batch_script( outdir, case, jobfile )
        elif ( cluster_type == "STAMPEDE" ): 
            write_and_submit_stampede2_slurm_script( outdir, case, jobfile )
        else: 
            queue_no = 1
            high_mem = True 
            write_and_submit_condor_script( outdir, case, executable, arguments, queue_no, high_mem )

def run_heterodimer_docking_calc( energy_fxn, rosetta_exe_path, cluster_type, test_set, restore, implicit_lipids, aqueous_pore ):
    """
    A function for running docking calculations on large heterodimeric sets

    Arguments
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = path to compiled Rosetta executable
        cluster_tupe = specify slurm or condor job submission
        test_set = Name of published set of test cases
        restore = Restore behavior to pre ref2015 for reference benchmark runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print("Initializing test: heterodimer docking")

    # Setup path to test sets and executables
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "mp_dock" + "." + platform + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking"
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + base_outdir )

    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    os.system( "mkdir " + test_outdir )

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Get chain info from test cases
        chains = case[5:7]
        pdbcode = case[0:4]
        partner_chains = chains[0] + "_" + chains[1]

        # Setup arguments by substitution
        native = inputs + "/" + test_set + "/" + case + "/" + case + ".pdb"
        prepacked = outdir + "/" + case + "_relaxed_prepacked.pdb"
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + ".span"
        scorefile = outdir + "/" + case + "_" + energy_fxn + "_docking.sc"
        s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners $partners -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 1000 -out:path:all $outdir -mp:dock:weights_fa $sfxn -out:file:scorefile $scfile" )
        arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir, partners=partner_chains, scfile=scorefile )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"
        if ( aqueous_pore == True ): 
            arguemnts = arguments + " -mp:pore:accomodate_pore true"

        # Write arguments and executable to a separate file
        jobfile = outdir + "/" + case + "_docking.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate job submission file and then submit to cluster
        print("Submitting docking test case from set " + test_set + ":", case)
        if ( cluster_type == "MARCC" ): 
            write_and_submit_slurm_batch_script( outdir, case, jobfile, str(10) )
        elif ( cluster_type == "STAMPEDE" ): 
            write_and_submit_stampede2_slurm_script( outdir, case, jobfile )
        else: 
            queue_no = 300
            write_and_submit_condor_script( outdir, case, executable, arguments, str(queue_no) )

def run_refined_native_heterodimer_docking_calc( energy_fxn, rosetta_exe_path, cluster_type, test_set, restore, implicit_lipids, aqueous_pore ):
    """
    A function for generating refined natives for docking calculations in the bound case (Hurwitz decoys)

    Arguments
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = path to compiled Rosetta executable
        cluster_tupe = specify slurm or condor job submission
        test_set = Name of published set of test cases
        restore = Restore behavior to pre ref2015 for reference benchmark runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print("Initializing test: heterodimer docking refined natives")

    # Setup path to test sets and executables
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "mp_dock" + "." + platform + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking"
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + base_outdir )

    # Generate a list of test cases from the input list file
    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Make a directory for the dataset of interest
    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    if ( not os.path.isdir( test_outdir ) ): 
        os.system( "mkdir " + test_outdir )

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        if ( not os.path.isdir( outdir ) ): 
            os.system( "mkdir " + outdir )

        # Make a directory for the refined natives
        refined_natives_dir = outdir + "/refined_natives"
        if ( not os.path.isdir( refined_natives_dir ) ): 
            os.system( "mkdir " + refined_natives_dir )

        # Change direcvtories into the location for refined natives
        os.system( "cd " + refined_natives_dir )

        # Get chain info from test cases
        chains = case[5:7]
        pdbcode = case[0:4]
        partner_chains = chains[0] + "_" + chains[1]

        # Setup arguments by substitution
        bound_native = inputs + "/" + test_set + "/" + case + "/" + case + ".pdb"
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + ".span"
        s = Template( " -in:file:s $native -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners $partners -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 10 -out:path:all $outdir -docking_local_refine -mp:dock:weights_fa $sfxn" )
        arguments = s.substitute( native=bound_native, spanfile=spanfile, sfxn=energy_fxn, outdir=refined_natives_dir, partners=partner_chains )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"
        if ( aqueous_pore == True ): 
            arguemnts = arguments + " -mp:pore:accomodate_pore true"

        # Write arguments and executable to a separate file
        jobfile = refined_natives_dir + "/" + case + "_refined_native_docking.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate job submission file and then submit to cluster
        print("Submitting docking test case from set " + test_set + ":", case)
        if ( cluster_type == "MARCC" ): 
            write_and_submit_slurm_batch_script( refined_natives_dir, case + "_refine_native", jobfile, str(10) )
        elif ( cluster_type == "STAMPEDE" ): 
            write_and_submit_stampede2_slurm_script( refined_natives_dir, case + "_refine_native", jobfile )
        else: 
            queue_no = 10
            write_and_submit_condor_script( refined_natives_dir, case + "_refine_native", executable, arguments, str(queue_no) )

def run_refined_native_docking_calc( energy_fxn, rosetta_exe_path, cluster_type, test_set, restore, implicit_lipids, aqueous_pore ): 
    """
    A function for generating refined natives for docking calculations in the bound case

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = path to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        test_set = Name of published set of test cases
        restore = Restore behavior to pre ref2015 for reference benchmark runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print("Initializing test: homodimer docking refined natives")

    # Setup path to test sets and executables
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "mp_dock" + "." + platform + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking"
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + base_outdir )

    # Generate a list of test cases from the input list file
    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Make a directory for the dataset of interest
    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + test_outdir )

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        if ( not os.path.isdir( outdir ) ): 
            os.system( "mkdir " + outdir )
        
        # Make a directory for the refined natives
        refined_natives_dir = outdir + "/refined_natives"
        if ( not os.path.isdir( refined_natives_dir ) ): 
            os.system( "mkdir " + refined_natives_dir )

        # Change directories into the location for refined natives
        os.system( "cd " + refined_natives_dir )

        # Setup arguments by substitution
        bound_native = inputs + "/" + test_set + "/" + case + "/" + case + "_AB_tr.pdb"
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + "_AB.span"
        s = Template( " -in:file:s $bound_native -in:file:native $bound_native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners A_B -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 10 -out:path:all $outdir -docking_local_refine -mp:dock:weights_fa $sfxn" )
        arguments = s.substitute( bound_native=bound_native, spanfile=spanfile, sfxn=energy_fxn, outdir=refined_natives_dir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"
        if ( aqueous_pore == True ): 
            arguemnts = arguments + " -mp:pore:accomodate_pore true"

        # Write arguments and executable to a separate file
        jobfile = refined_natives_dir + "/" + case + "_make_bound_refined.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate job submission file and then submit to cluster
        print("Submitting bound refine native job for case from set " + test_set + ":", case)
        if ( cluster_type == "MARCC" ): 
            write_and_submit_slurm_batch_script( refined_natives_dir, case + "_refined_native", jobfile, str(10) )
        elif ( cluster_type == "STAMPEDE" ): 
            write_and_submit_stampede2_slurm_script( refined_natives_dir, case + "_refined_native", jobfile )
        else: 
            queue_no = 10
            write_and_submit_condor_script( refined_natives_dir, case + "_refined_native", executable, arguments, str(queue_no) )
 
def run_docking_calc( energy_fxn, rosetta_exe_path, cluster_type, test_set, restore, implicit_lipids, aqueous_pore ): 
    """
    A function for running the docking calculations needed for the docking test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = path to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        test_set = Name of published set of test cases
        restore = Restore behavior to pre ref2015 for reference benchmark runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print("Initializing test: homodimer docking")

    # Setup path to test sets and executables
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "mp_dock" + "." + platform + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking"
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + base_outdir )

    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    os.system( "mkdir " + test_outdir )

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup arguments by substitution
        native = inputs + "/" + test_set + "/" + case + "/" + case + "_AB_tr.pdb"
        prepacked = inputs + "/" + test_set + "/" + case + "/" + case + "_AB_tr.prepack.pdb"
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + "_AB.span"
        s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners A_B -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct 1000 -out:path:all $outdir" )
        arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"
        if ( aqueous_pore == True ): 
            arguemnts = arguments + " -mp:pore:accomodate_pore true"

        # Write arguments and executable to a separate file
        jobfile = outdir + "/" + case + "_seqrecov.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate job submission file and then submit to cluster
        print("Submitting docking test case from set " + test_set + ":", case)
        if ( cluster_type == "MARCC" ): 
            write_and_submit_slurm_batch_script( outdir, case, jobfile, str(10) )
        elif ( cluster_type == "STAMPEDE" ): 
            write_and_submit_stampede2_slurm_script( outdir, case, jobfile )
        else: 
            queue_no = 300
            write_and_submit_condor_script( outdir, case, executable, arguments, str(queue_no) )

##@brief: Stage 1: Generate 10 refined natives 
def make_relaxed_homology_models( energy_fxn, test_set, rosetta_exe_path, restore, implicit_lipids ): 
    """
	A function for refining homology models in their native context

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        restore = Option to restore talaris behavior before ref2015 for reference runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print( "Creating refined homology models in the resident energy function" )

    # Setup paths to input files and test executable
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "rosetta_scripts." + platform + compiler + buildenv
    xml_script = benchmark + "xml/test_3.3_decoy_refinement.xml"

    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Create paths to output files
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/"
    os.system( "mkdir " + base_outdir )
    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    os.system( "mkdir " + test_outdir )

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Make a directory for the refined models
        models_outdir = outdir + "/refined_models"
        os.system( "mkdir " + models_outdir ) 

        # Setup arguments by substitution
        native = inputs + "/" + test_set + "/" + case + "/" + case + ".pdb"
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + ".span"
        scorefile = case + "_" + energy_fxn + "_refined_models.sc"
        s = Template( "-in:file:native $native -in:file:s $native -mp:setup:spanfiles $spanfile -parser:script_vars sfxn_weights=$energy_fxn -parser:protocol $xml -out:file:scorefile $scorefile -out:path:all $outdir -nstruct 10 -run:multiple_processes_writing_to_one_directory")
        arguments = s.substitute( native=native, spanfile=spanfile, scorefile=scorefile, xml=xml_script, energy_fxn=energy_fxn, outdir=models_outdir )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"

        # Write arguments and executable to a separate file
        jobfile = models_outdir + "/" + case + "_refine_models.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )

        # Generate job submission file and then submit to cluster
        print("Submitting docking test case from set " + test_set + ":", case)
        queue_no = 10
        write_and_submit_condor_script( models_outdir, case, executable, arguments, str(queue_no) )

##@brief: Stage 2 - Identify the lowest scoring relaxed case and make a copy
def find_lowest_scoring_relaxed_native( energy_fxn, test_set ): 
    """
	A function for finding the lowest scoring relaxed model

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
		test_set = Name of test set containing target models
    """

    # Setup paths to input files and test executable
    inputs = benchmark + "inputs/test_3.2_docking"

	# Establish path to test cases
    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set
    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
    	test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

		# Navigate to the directory contianing output files
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
        models_dir = outdir + "/refined_models"
        os.chdir( models_dir )

		# Remove the first line of the scorefile
        scorefile_orig = case + "_" + energy_fxn + "_refined_models.sc"
        scorefile_no_seq = case + "_" + energy_fxn + "_refined_models_no_seq.sc"
        os.system( "tail -n +3 " + scorefile_orig + " > " + scorefile_no_seq )

		# Sort the file contents by total score
        scorefile_sorted = case + "_" + energy_fxn + "_sorted.sc"
        os.system( "sort -k2 " + scorefile_no_seq + " > " + scorefile_sorted )

		# Read in contents of the file and extract the description of the lowest scoring
        with open( scorefile_sorted, 'r' ) as f: 
            contents = f.readlines()
            contents = [ x.strip() for x in contents ]
            contents = [ x.split() for x in contents ]

        lowest_desc = models_dir + "/" + contents[-1][-1] + ".pdb" 
        lowest_copy = models_dir + "/" + case + "_relaxed_lowest.pdb"

        # Make a copy of this file to mark it as the lowest
        os.system( "cp " + lowest_desc + " " + lowest_copy )

##@brief: Stage 3: Generate a prepacked structure for the lowest scoring native
def make_prepacked_docking_input( energy_fxn, test_set, rosetta_exe_path, restore, implicit_lipids ): 
    """
	A function for creating a prepacked structure for docking

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        restore = Option to restore talaris behavior before ref2015 for reference runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print( "Generating a prepacked input structure from the lowest scoring relaxed model" )

    # Setup path to inputs and executables
    inputs = benchmark + "inputs/test_3.2_docking"
    executable = rosetta_exe_path + "docking_prepack_protocol" + "." + platform + compiler + buildenv
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking"
    if ( not os.path.isdir( base_outdir ) ): 
        os.system( "mkdir " + base_outdir )

    # Establish path to output models
    test_outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set

    # Curate a list of test cases
    list_of_test_cases = inputs + "/" + test_set + "/dimers.list" 
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # For each test case, generate specific arguments, job files, and then run
    for case in test_cases:

        # Make one directory per case
        outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/" + test_set + "/" + case
       	if ( not os.path.isdir( outdir ) ):
	        os.system( "mkdir " + outdir )
        os.chdir( outdir )

        # Move the lowest scoring relaxed model to the main directory
        lowest_scoring = case + "_relaxed_lowest.pdb"
        os.system( "cp refined_models/" + lowest_scoring + " " + lowest_scoring )

        # Get chain info from test cases
        partner_chains = "A_B"
        if ( len(case) > 4 ): 
	        chains = case[5:7]
	        pdbcode = case[0:4]
	        partner_chains = chains[0] + "_" + chains[1]

        # Setup arguments by substitution
        native = lowest_scoring
        spanfile = inputs + "/" + test_set + "/" + case + "/" + case + ".span"
        scorefile = case + "_prepacking_step.sc"
        s = Template( " -overwrite -in:file:s $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -docking:partners $partners -packing:pack_missing_sidechains 0 -nstruct 1 -out:path:all $outdir -out:file:scorefile $scfile" )
        arguments = s.substitute( native=native, spanfile=spanfile, sfxn=energy_fxn, outdir=outdir, partners=partner_chains, scfile=scorefile )
        if ( restore == True ): 
            arguments = arguments + " -restore_talaris_behavior"
        if ( implicit_lipids == True ): 
            arguments = arguments + " -mp:lipids:use_implicit_lipids true -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"

        # Write arguments and executable to a separate file
        jobfile = outdir + "/" + case + "_prepack_docking_model.sh"
        with open( jobfile, 'a' ) as f: 
            f.write( "#!/bin/bash\n" )
            f.write( executable + " " + arguments + "\n" )
            f.close()
        os.system( "chmod +x " + jobfile )
        os.system( "bash " + jobfile )

        # Rename the resulting PDB file as input
        os.system( "cp " + case + "_relaxed_lowest_0001.pdb " + case + "_relaxed_prepacked.pdb" )
        
def run_decoy_discrimination_calc( energy_fxn, rosetta_exe_path, cluster_type, restore, implicit_lipids, aqueous_pore ): 
    """
    A function for running the docking calculations needed for the docking test

    Arguments: 
        energy_fxn = energy function to use (typically, name of the weights file)
        rosetta_exe_path = pah to compiled Rosetta executable
        cluster_type = specify slurm or condor job submission
        restore = Option to restore talaris behavior before ref2015 for reference runs
        implicit_lipids = Use implicit lipid parameters? 
        aqueous_pore = Account for an aqueous pore? 
    """

    print("Initializing test: decoy-discrimination")

    ### Setup some general variables
    inputs = benchmark + "inputs/test_3.3_decoy_discrimination"
    executable = rosetta_exe_path + "rosetta_scripts." + platform + compiler + buildenv
    xml_script = benchmark + "xml/test_3.3_decoy_refinement.xml"

    ### Make the base decoy discrimination output directory
    base_outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination"
    os.system( "mkdir " + base_outdir )

    ### Test Set #1: Yarov-Yaravoy Low Resolution Decoys (Membrane ab initio generated)
    #list_of_test01_cases = inputs + "/yarov-yaravoy-set/decoy_sets.list"
    #with open( list_of_test01_cases, 'rb' ) as f:
    #    test01_cases = f.readlines()
    #test01_cases = [ x.strip() for x in test01_cases ]

    #outdir_test01 = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/yarov-yaravoy-set"
    #os.system( "mkdir " + outdir_test01 )
    #os.chdir( outdir_test01 )

    # For each test case, generate specific arguments, condor_files, and then run
    #for case in test01_cases:

    #    outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/yarov-yaravoy-set/" + case
#        os.system( "mkdir " + outdir )
#        os.system( "cd " + outdir )

        # Select 100 random decoys of the set of 5000
        #decoy_selection = random.sample(range(1,5000), 100)
        #decoy_selection = [ str(x).zfill(4) for x in decoy_selection ]

        # Open a file to determine the prefix
        #prefixfile = inputs + "/yarov-yaravoy-set/" + case + "/prefix.txt"
        #with open( prefixfile, 'rb' ) as f: 
        #    contents = f.readlines()
        #prefix = contents[0].strip()

        # Make a list of PDBs to refine
        #pdblist = outdir + "/random_models.list"
        #with open( pdblist, 'w' ) as f: 
        #    for decoy in decoy_selection: 
        #        pdb = inputs + "/yarov-yaravoy-set/" + case + "/" + prefix + decoy + ".pdb"
         #       f.write( pdb + "\n" )
         #   f.close()

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
 #       native = inputs + "/yarov-yaravoy-set/" + case + "/" + case + "_native.pdb"
 #       spanfile = inputs + "/yarov-yaravoy-set/" + case + "/" + case + ".span"
 #       pdblist = inputs + "/yarov-yaravoy-set/" + case + "/" + case + "_sampled_candidates.list"
 #       s = Template( "-overwrite -in:file:native $native -relax:constrain_relax_to_start_coords -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
 #       arguments = s.substitute( modellist=pdblist, span=spanfile, xml=xml_script, sfxn=Options.energy_fxn, native=native, outdir=outdir)
 #       if ( restore == True ): 
 #           arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
 #       if ( implicit_lipids == True ): 
 #           arguments = arguments + " -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"
    
        # Write arguments and executable to a separate file
 #       jobfile = outdir + "/" + case + "_seqrecov.sh"
 #       with open( jobfile, 'a' ) as f: 
 #           f.write( "#!/bin/bash\n" )
 #           f.write( executable + " " + arguments + "\n" )
 #           f.close()
 #       os.system( "chmod +x " + jobfile )

        # Generate a condor submission file and submit the job to Jazz
  #      condor_case_name = case + "_decoy_disc"
  #      print("Submitting decoy-discrimination test case from Yarov-Yaravoy set:", condor_case_name)
  #      if ( cluster_type == "MARCC" ): 
  #          write_and_submit_slurm_batch_script( outdir, condor_case_name, jobfile )
  #      elif ( cluster_type == "STAMPEDE" ): 
  #          write_and_submit_stampede2_slurm_script( outdir, condor_case_name, jobfile )
  #      else: 
  #          write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )


    ### Test Set #2: Dutagaci High Resolution Decoys (Molecular dynamics generated)
    list_of_test02_cases = inputs + "/dutagaci-set/decoy_sets.list"
    with open( list_of_test02_cases, 'rb' ) as f:
        test02_cases = f.readlines()
    test02_cases = [ x.strip() for x in test02_cases ]

    outdir_test02 = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/dutagaci-set"
    os.system( "mkdir " + outdir_test02 )
    os.chdir( outdir_test02 )

    # For each test case, generate specific arguments, condor_files, and then run
    for case in test02_cases:

        # Open list of decoy lists
        with open( inputs + "/dutagaci-set/" + case + "/list.of.decoy.lists", 'rt' ) as f: 
            list_of_decoy_lists = f.readlines()
            list_of_decoy_lists = [ x.strip() for x in list_of_decoy_lists ] 

        outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination/dutagaci-set/" + case
        os.system( "mkdir " + outdir )
        os.system( "cd " + outdir )

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        native = inputs + "/dutagaci-set/" + case + "/" + case + "_native.pdb"
        spanfile = inputs + "/dutagaci-set/" + case + "/" + case + ".span"
        
        # Iterate over all decoy lists
        i = 1
        for decoy_set in list_of_decoy_lists: 

            modelslist = inputs + "/dutagaci-set/" + case + "/" + decoy_set
            s = Template( " -overwrite -in:file:native $native -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
            arguments = s.substitute( modellist=modelslist, span=spanfile, xml=xml_script, sfxn=Options.energy_fxn, native=native, outdir=outdir)
            if ( restore == True ): 
                arguments = arguments + " -restore_talaris_behavior -restore_lazaridis_imm_behavior"
            if ( implicit_lipids == True ): 
                arguments = arguments + " -mp:lipids:temperature 37.0 -mp:lipids:composition DLPC"

            # Write arguments and executable to a separate file
            jobfile = outdir + "/" + case + "_relax_batch_" + str(i) + ".sh"
            with open( jobfile, 'a' ) as f: 
                f.write( "#!/bin/bash\n" )
                f.write( executable + " " + arguments + "\n" )
                f.close()
            os.system( "chmod +x " + jobfile )

            # Generate a condor submission file and submit the job to Jazz
            condor_case_name = case + "_dutagaci_" + str(i)
            print("Submitting decoy-discrimination test case from Dutagaci set:", condor_case_name)
            if ( cluster_type == "MARCC" ): 
                write_and_submit_slurm_batch_script( outdir, condor_case_name, jobfile )
            elif ( cluster_type == "STAMPEDE" ): 
                write_and_submit_stampede2_slurm_script( outdir, condor_case_name, jobfile )
            else: 
                write_and_submit_condor_script( outdir, condor_case_name, executable, arguments )

            i = i+1

def main( args ):

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_t01 --which_tests all" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function weights file", )

    parser.add_option('--stable', '-s', 
        action="store", 
        help="For reference runs, use the stable Rosetta master branch", )

    parser.add_option('--cluster_type', '-t', 
        action="store", 
        help="Specify option to run jobs on a slurm or condor system", )

    parser.add_option( '--which_tests', '-w', 
        action="store", 
        help="Specify which test groups to run. Options are: ddG, landscape, prediction", )

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
    if ( not Options.energy_fxn or not Options.stable or not Options.cluster_type ): 
        print("Missing required options --energy_fxn, --stable, and/or --cluster_type! Exiting...")
        sys.exit()

    if ( ( not Options.cluster_type == "MARCC" ) and ( not Options.cluster_type == "CONDOR" ) and ( not Options.cluster_type == "STAMPEDE") ): 
        print("Invalid option for --cluster_type. Currently only supporting STAMPEDE, MARCC, or CONDOR")
        sys.exit()

    test_types = []
    if ( not Options.which_tests ): 
        test_types.append( "ddG" )
        test_types.append( "landscape" )
        test_types.append( "prediction" )
    else: 
        which_types = Options.which_tests.split(",")
        for t in which_types: 
            if ( t == "ddG" or t == "landscape" or t == "prediction" ): 
                test_types.append( t )
            else: 
                print("Invalid test type", t, "Exiting...")
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

    # Run energy landscape testing group
    if ( "landscape" in test_types ): 
    
        # Energy landscape test for single TM peptides found in nature
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.1_monomer_landscape", "helices.list", "xml/test_1.1_monomer_landscape.xml", restore, include_lipids )

        # Energy landscape test for aromatic-capped peptides
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.2_aro_landscape", "aro_helices.list", "xml/test_1.2_aro_landscape.xml", restore, include_lipids, "true" )

        # Energy landscape test for leucine-lysine peptides
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.3_lk_landscape", "lk_peptides.list", "xml/test_1.3_lk_landscape.xml", restore, include_lipids, "true" )

        # Energy landscape test for adsorbed peptides
        run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_1.4_adsorbed_pept_landscape", "adsorbed_peptides.list", "xml/test_1.3_lk_landscape.xml", restore, include_lipids, "true" )

    # Run prediction calculations
    if ( "prediction" in test_types ): 

        # Fixed backbone design calculation for sequence recovery test
        #run_fixed_backbone_design_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, restore, include_lipids, add_pore )

        # Decoy Discrimination calculation for large and small sets
        run_decoy_discrimination_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, restore, include_lipids, False )

        # Docking calculation for small homodimer set (Lomize et al. 2017)
        #run_refined_native_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "small-homodimer-set", restore, include_lipids, False )
        #run_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "small-homodimer-set", restore, include_lipids, False )
        #make_relaxed_homology_models( Options.energy_fxn, "small-homodimer-set", rosetta_exe_path, restore, include_lipids )

        # Docking calculation for large homodimer set (Alford & Koehler Leman 2015)
        #run_refined_native_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-homodimer-set", restore, include_lipids, False )
        #run_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-homodimer-set", restore, include_lipids, False )
        #make_relaxed_homology_models( Options.energy_fxn, "large-homodimer-set", rosetta_exe_path, restore, include_lipids )

        # Docking calculation for large bound-bound set (Hurwitz et al. 2016)
        #run_refined_native_heterodimer_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-bound-set", restore, include_lipids, False )
        #run_heterodimer_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-bound-set", restore, include_lipids, False )
        #make_relaxed_homology_models( Options.energy_fxn, "large-bound-set", rosetta_exe_path, restore, include_lipids )

        # Docking calculation for large unbound set (simulated, Hurwitz et al. 2016)
        #make_relaxed_homology_models( Options.energy_fxn, "large-unbound-set", rosetta_exe_path, restore, include_lipids )
        #find_lowest_scoring_relaxed_native( Options.energy_fxn, "large-unbound-set" )
        #make_prepacked_docking_input( Options.energy_fxn, "large-unbound-set", rosetta_exe_path, restore, include_lipids )
        #run_heterodimer_docking_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "large-unbound-set", restore, include_lipids, False )
    
    # Run ddG calculations
    if ( "ddG" in test_types ): 

        # ddG of insertion landscape calculation for Ulmschneider set
        #run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_2.2_ddG_of_insertion", "insertion_peptide.dat", "xml/test_2.2_ddG_insertion_landscape.xml", restore, include_lipids, add_pore, "true" )

        # ddG of insertion landscape calculation for pH dependent set - generate at pH = 4
        #run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_2.3_pH_dependent_insertion", "pH-inserted-helices.list", "xml/test_2.3_pH_landscape.xml", restore, include_lipids, add_pore, "true", "4" )
    
        # ddG of insertion landscape calculation for pH dependent set - generate at pH = 7
        #run_energy_landscape_calc( Options.energy_fxn, rosetta_exe_path, Options.cluster_type, "test_2.3_pH_dependent_insertion", "pH-inserted-helices.list", "xml/test_2.3_pH_landscape.xml", restore, include_lipids, add_pore, "true", "7" )

        # ddG of mutation calculation for Moon & Fleming Set
        #run_ddG_of_mutation_calc( Options.energy_fxn, "OmpLA/OmpLA_Moon_Fleming_set.dat", "OmpLA_Moon_Fleming_set", restore,  include_lipids, add_pore )

        # ddG of mutation calculation for McDonald & Fleming Set
        #run_ddG_of_mutation_calc( Options.energy_fxn, "OmpLA_aro/OmpLA_aro_McDonald_Fleming_set.dat", "OmpLA_aro_McDonald_Fleming_set", restore, include_lipids, add_pore )

        # ddG of mutation calculation for Marx & Fleming set
        run_ddG_of_mutation_calc( Options.energy_fxn, "PagP/PagP_Marx_Fleming_set.dat", "PagP_Marx_Fleming_set", restore, include_lipids, False )

if __name__ == "__main__" : main(sys.argv)
