#!/usr/bin/env python
###################################################################
#@file:         make_protein_energy_landscape.py                                                                                   
#@description:  Generate orientation-dependent peptide energy map                                               
#@args:			    --energy_fxn (wts) --targets (list of targets)                                                          
#@author: 		  Rebecca F. Alford                   
#@email: 		    rfalford12@gmail.com                                          
###################################################################

import random
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_protein_energy_landscape_calc( energy_fxn, config, targets, test_name, xml ): 
	"""
	A function to generate orientation-dependent peptide energy landscapes

	Arguments: 
		energy_fxn = energy function to use for calculations (typically, name of the weights file)
		config = path to benchmark, Rosetta executables
		targets = list of targets
		xml = Path to RosettaScript defining the peptide landscape search protocol
	"""

    print( "Making energy landscapes for the multipass membrane protein set", targets ) 

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

    # For each test case, generate specific arguments, a condor file, and then run
    for case in test_cases:

        # Setup case-specific variables (pdbfile, spanfile, xmlargs)
        pdbfile = targets + "/" + case + "/" + case + ".pdb"
        spanfile = targets + "/" + case + "/" + case + ".span"

        # Default composition is DLPC
        # Should I tune the pH of my simulation?
        arguments = " -overwrite -in:file:s " +  pdbfile + " -mp:setup:spanfiles " + spanfile + " -parser:script_vars sfxn_weights=" + energy_fxn + " -parser:protocol " + xml_script 
        arguments = arguments + " -mp:lipids:composition DLPC "

  	casedir = outdir + "/" + case
    if ( not os.path.isdir( casedir ) ): 
      os.system( "mkdir " + casedir )
      os.chdir( casedir )

    # Write jobfile and submit to the HPC
    print("Submitting orientation sampling calculations for energy landscape case:", case ) 
    jobname = case + "_energy_landscape"
    jobfile = hpc_util.make_jobfile( casedir, case, executable, arguments )
    hpc_util.submit_condor_job( config.benchmark_path, jobname, jobfile, "" )
