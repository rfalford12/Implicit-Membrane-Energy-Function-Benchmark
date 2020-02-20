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

import sys, os, random, read_config, hpc_util
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

def convert_spanfile_to_string( spanfile ): 

	with open( spanfile, 'rt' ) as f: 
		spans = f.readlines()
		spans = [ x.strip() for x in spans ]

	# Parse tm regions, skipping the first four header lines
	span_info = []
	for tm_span in range(4, len(spans)): 
		single_span_data = spans[tm_span].split( ' ' )
		single_span_info = []
		for entry in single_span_data: 
			if ( entry != "" ): 
				single_span_info.append( entry )
		span_info.append( single_span_info )

	# Convert list of spans into a string that will be used as input for kink finder
	spanstr = ""
	for span in span_info: 
		spanstr = spanstr + span[0] + "-" + span[1] + " "

	return spanstr

def main( args ): 

	# Read options from the command line
	parser = OptionParser(usage="usage %prog --energy_fxn franklin2019 --which_tests all" )
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


	# Test #8: Sequence recovery calculation
	if ( "sequence-recovery" in test_names ): 

		# Make list of native and designed PDB files
		datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/sequence-recovery/"
		os.chdir( datadir )
		os.system( "ls */*_0001.pdb > designed.list" )
		with open( "designed.list",'rt' ) as f: 
			contents = f.readlines()
			contents = [ x.strip() for x in contents ]
			pdbid = [ x.split("/")[0] for x in contents ]

		with open( "natives.list", 'wt' ) as f: 
			basedir = config.benchmark_path + "targets/design/"
			for pdb in pdbid: 
				pdbpath = basedir + pdb + "/" + pdb + "_tr_ignorechain.pdb\n"
				f.write( pdbpath )

		# Run mp_seqrecov application
		executable = config.rosetta_path + "mp_seqrecov." + config.platform + config.compiler + config.buildenv
		output_file = Options.energy_fxn + "_seqrecov.txt"
		s = Template( " -overwrite -native_pdb_list $natives -redesign_pdb_list $designed -seq_recov_filename $outfile -in:ignore_unrecognized_res -read_only_ATOM_entries")
		arguments = s.substitute( natives="natives.list", designed="designed.list", outfile=output_file )
		if ( restore == False ): 
			arguments = arguments + " -mp:lipids:composition DLPC -mp:lipids:temperature 37"
		os.system( executable + arguments )

		# Process the sequence recovery data
		os.system( "python tests/python/process_protein_design_results.py --energy_fxn " + Options.energy_fxn + " --basedir " + config.benchmark_path + " --seqrecov_file " + output_file )

	# Test #9: Side chain distribution calculations
	if ( "sc-distribution" in test_names ): 

		# Check for existence of designed and native lists
		datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/sequence-recovery/"
		os.chdir( datadir )
		if ( not os.path.isfile( "natives.list") or not os.path.isfile( "designed.list" ) ): 
			os.system( "ls */*_0001.pdb > designed.list" )
			with open( "designed.list",'rt' ) as f: 
				contents = f.readlines()
				contents = [ x.strip() for x in contents ]
				pdbid = [ x.split("/")[0] for x in contents ]

			with open( "natives.list", 'wt' ) as f: 
				basedir = config.benchmark_path + "targets/design/"
				for pdb in pdbid: 
					pdbpath = basedir + pdb + "/" + pdb + "_tr_ignorechain.pdb\n"
					f.write( pdbpath )

		# Run predict side chain distribution script
		predict_side_chain_distribution.compute_side_chain_distribution( config, datadir + "natives.list", datadir + "/designed.list" )

	if ( "decoy-discrimination" in test_names ): 

		# for each target, make a list of pdb decoys
		targets = [ "brd7", "fmr5", "ltpa", "rhod", "vatp"]
		datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/decoy-discrimination/"
		# first hires targets
		hires_dir = datadir + "hires"
		os.chdir( hires_dir )
		for target in targets: 
			os.chdir( target )
			os.system( "ls decoy*.pdb > decoys.list" )
			os.chdir( "../" )

		# then lowres targets
		lowres_dir = datadir + "lowres"
		os.chdir( lowres_dir )
		for target in targets: 
			os.chdir( target )
			os.system( "ls *.pdb > decoys.list" )
			os.chdir( "../" )

		# rescore each to calculate the rms and total score
		executable = config.rosetta_path + "score_jd2" + "." + config.platform + config.compiler + config.buildenv
		s = Template( "-in:file:l decoys.list -in:file:native $native -mp:setup:spanfiles from_structure -out:file:scorefile $scorefile -in:membrane")
		for target in targets: 

			# hires
			targetdir = hires_dir + "/" + target
			os.chdir( targetdir )
			output_scores = target + "_hires.sc"
			native_pdb = config.benchmark_path + "targets/structure/D6_decoy_discrimination/hires/" + target + "/" + target + "_native.pdb"
			spanfile = config.benchmark_path + "targets/structure/D6_decoy_discrimination/hires/" + target + "/" + target + ".span"
			arguments = s.substitute( scorefile=output_scores, native=native_pdb )
			jobname = "rescore_hires_" + target
			hpc_util.submit_condor_job( targetdir, jobname, executable, arguments, 1 )

			# lowres
			targetdir = lowres_dir + "/" + target
			os.chdir( targetdir )
			output_scores = target + "_lowres.sc"
			native_pdb = config.benchmark_path + "targets/structure/D6_decoy_discrimination/lowres/" + target + "/" + target + "_native.pdb"
			spanfile = config.benchmark_path + "targets/structure/D6_decoy_discrimination/lowres/" + target + "/" + target + ".span"
			arguments = s.substitute( scorefile=output_scores, native=native_pdb )
			jobname = "rescore_lowres_" + target
			hpc_util.submit_condor_job( targetdir, jobname, executable, arguments, 1 )

	if ( "helix-kinks" in test_names ): 

		# Read list of test case IDs and PDBs
		targets_path = config.benchmark_path + "targets/structure/D5_helix_kinks/"
		list_of_targets = targets_path + "targets.list"
		with open( list_of_targets, 'rt' ) as f: 
			targets = f.readlines()
			targets = [ x.strip() for x in targets ]

		# Change directories to a data analysis dir
		outdir = config.benchmark_path + "data/" + Options.energy_fxn + "/helix-kinks/" 
		if ( not os.path.isdir( outdir ) ): 
			aya.exit( "Data generation for helix kinks was skipped! Cannot process non-existent data. Exiting!" )
		os.chdir( outdir )

		# Iterate through each target
		for target in targets:

			# Clean up the directories with refined models
			casedir = outdir + target 
			os.chdir( casedir )
			os.system( "rm *.in_progress" )
			os.system( "mkdir metadata" )
			os.system( "mv *.out metadata/" )
			os.system( "mv *.err metadata/" )
			os.system( "mv *.log metadata/" )
			os.system( "mkdir scorefiles" )
			os.system( "mv *.sc scorefiles/" )

			# Make a list of the refined models
			os.system( "ls *.pdb > models.list" )
			with open( "models.list" ) as f: 
				refined_models = f.readlines()
				refined_models = [ x.strip() for x in refined_models ]

			# Conver tthe spanfile from the targets directory to a string format
			kink_finder_script = config.benchmark_path + "external/kink-finder/Kink_Finder.py"
			spanfile = config.benchmark_path + "targets/structure/D5_helix_kinks/" + target + "/" + target + ".span"
			spanstr = convert_spanfile_to_string( spanfile )
			os.system( "mkdir kink_data" )
			i = 1
			#for m in (refined_models): 
			kink_cmd = "python2 " + kink_finder_script + " -f " + refined_models[0] + " -l '" + spanstr + "'"
			os.system( kink_cmd )
			os.system( "mv kinks.csv kink_data/" + target + "_" + str(i) + ".csv" ) 

			# Run kink finder on each model
			# will need to re-format the helix definitions here and then run kink finder
			# on each model, output into a kink files directory

			# then run the process kink data script

	if ( "protein-protein-docking" in test_names ): 

		print("temp")
		# Analyze docking decoys
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", False, True )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", False  )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", False  )

		# Analyze local refine decoys
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D2_single_TM_complexes", "protein-protein-docking", True, True )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D3_multi_TM_bound_complexes", "protein-protein-docking", True  )
		make_asymm_docked_complexes.analyze_interfaces( Options.energy_fxn, config, "D4_multi_TM_unbound_complexes", "protein-protein-docking", True  )


if __name__ == "__main__" : main(sys.argv)

