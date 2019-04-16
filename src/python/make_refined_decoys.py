#!/usr/bin/env python
###################################################################
#@file:         make_refined_models.py                                  
#@description:  Generate relaxed candidate structures
#@args:       --energy_fxn (wts) --targets (list of targets)                                                          
#@author:     Rebecca F. Alford                   
#@email:    rfalford12@gmail.com                                          
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def run_refinement_calc( energy_fxn, config, targets, test_name, xml, nstruct, subset="" ): 
  """
    A function for refining canddiate structures for the decoy discrimination test

    Arguments: 
      energy_fxn = energy function to use for calculations (typically, name of the weights file)
      config = path to benchmark, Rosetta executables
      targets = list of targets
      test_name = Name of test
  """

  print "Refining candidate structures for dataset", targets

  # Read list of test case IDs and PDBs
  list_of_targets = config.benchmark_path + "targets/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
      temp = f.readlines()
      temp = [ x.strip() for x in temp ]
      target_ids = [ x.split("/")[0] for x in temp ]
      target_pdbs = [ x.split("/")[1] for x in temp]

  # Generate path to executable
  executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
  xml_script =  config.benchmark_path + "src/xml/" + xml

  # Change directories to a data analysis dir
  outdir = config.benchmark_path + "data/" + energy_fxn + "/" + test_name 
  if ( not os.path.isdir( outdir ) ): 
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

  if ( subset != "" ): 
    outdir = outdir + "/" + subset
    if ( not os.path.isdir( outdir ) ): 
      os.system( "mkdir " + outdir )
      os.chdir( outdir )

  for i in range(0, len(target_ids)): 

    # Setup case-specific variables (pdbfile, spanfile)
    native = targets + "/" + target_ids[i] + "/" + target_pdbs[i]
    spanfile = targets + "/" + target_ids[i] + "/" + target_ids[i] + ".span"

    casedir = outdir + "/" + target_ids[i]
    if ( not os.path.isdir( casedir ) ): 
      os.system( "mkdir " + casedir )
      os.chdir( casedir )

    # Generate a string of arguments from the case-specific variables
    s = Template( "-overwrite -in:file:native $native -relax:constrain_relax_to_start_coords -in:file:s $model -mp:setup:spanfiles $span -parser:script_vars sfxn_weights=$sfxn -parser:protocol $xml -out:file:scorefile refined_models.sc -out:path:all $outdir -nstruct $nmodels")
    arguments = s.substitute( model=native, span=spanfile, xml=xml_script, sfxn=energy_fxn, native=native, outdir=casedir, nmodels=nstruct)

    # Write jobfile and submit to the HPC
    print "Submitting refinement calculations for case:", target_ids[i] 
    jobname = target_ids[i] + "_refine"
    hpc_util.submit_condor_job( casedir, jobname, executable, arguments )
