#!/usr/bin/env python
###################################################################
#@file:         make_symm_docked_complexes.py                                                                                   
#@description:  Assemble symmetric membrane protein complexes                                             
#@args:			    --energy_fxn (wts) --targets (list of targets)                                                          
#@author: 		  Rebecca F. Alford                   
#@email: 		    rfalford12@gmail.com                                          
###################################################################

import random
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_symmdock_calc( energy_fxn, config, targets, test_name ): 
	"""
		A function for assembling symmetric membrane protein complexes

		Arguments: 
			energy_fxn = energy function to use for calculations (typically, name of the weights file)
			config = path to benchmark, Rosetta executables
			targets = list of targets
			test_name = Name of test
	"""

	print( "Refining candidate structures for dataset", targets ): 

   	# Read list of energy landscape test cases
    list_of_targets = config.benchmark_path + "targets/" + targets + "/targets.list"
    with open( list_of_targets, 'rt' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    # Generate path to executable
    executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
    xml_script =  config.benchmark_path + "/src/xml/" + xml_protocol

    # Change directories to a data analysis dir
    outdir = config.benchmark_path + "/data/" + energy_fxn + "/" + test_name 
    if ( not os.path.isdir( outdir ) ): 
    	os.system( "mkdir " + outdir )
    	os.chdir( outdir )

    for case in test_cases: 

    	# Setup case-specific variables (pdbfile, spanfile)
    	native = targets + "/" + case + "/" + case + "_native.pdb"
    	spanfile = targets + "/" + case + "/" + case + ".span"
    	pdblist = targets + "/" + case + "/" + case + ".decoys.list"

    	casedir = outdir + "/" + case
  		if ( not os.path.isdir( casedir ) ): 
  			os.system( "mkdir " + casedir )
  			os.chdir( casedir )

    	# Generate a string of arguments from the case-specific variables
		s = Template( "-overwrite -in:file:native $native -relax:constrain_relax_to_start_coords -in:file:l $modellist -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir")
 		arguments = s.substitute( modellist=pdblist, span=spanfile, xml=xml_script, sfxn=energy_fxn, native=native, outdir=casedir)

		# Write jobfile and submit to the HPC
  		print("Submitting refinement calculations for decoy disc case:", case ) 
  		jobname = case + "_refine"
		jobfile = hpc_util.make_jobfile( casedir, case, executable, arguments )
		hpc_util.submit_condor_job( config.benchmark_path, jobname, executable, arguments )
