#!/usr/bin/env python
###################################################################
#@file:         make_asymm_docked_complexes.py                                                                                   
#@description:  Assemble membrane protein complexes                                             
#@args:         --energy_fxn (wts) --targets (list of targets)                                                          
#@author:       Rebecca F. Alford                   
#@email:        rfalford12@gmail.com                                          
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )


def post_process_prepack_pdb( energy_fxn, config, targets, test_name, transformed=False ): 
  """
    A function for generating prepacked structures for membrane protein-protein docking

    Arguments: 
      energy_fxn = energy function to use for calculations (typically, name of the weights file)
      config = path to benchmark, Rosetta executables
      targets = list of targets to use (homodimers, hetero-dimers unbound, hetero-dimers bound)
    """

  print("Generating prepacked structures for protein-protein complex set", targets)

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/structure/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
      temp = f.readlines()
      temp = [ x.strip() for x in temp ]
      target_ids = [ x.split("/")[0] for x in temp ]
      target_pdbs = [ x.split("/")[1] for x in temp ]

  # Generate path to executable
  executable = config.rosetta_path + "docking_prepack_protocol" + "." + config.platform + config.compiler + config.buildenv

  for i in range(0, len(target_ids)): 

    # Change directories to a data analysis dir
    outdir = config.benchmark_path + "/data/" + energy_fxn + " /protein-protein-docking/" + targets + "/" + target_ids[i]
    if ( not os.path.isdir( outdir ) ): 
      sys.exit( "Prepacking step was skipped! Cannot clean non-existent prepack structures" )
    os.chdir( outdir )

    # Setup case-specific variables (pdbfile, spanfile)
    partners = target_pdbs[i].split("_")[1][0] + target_pdbs[i].split("_")[1][1]
    prepack_pdb_title = target_ids[i] + "_0001"
    if ( transformed ):  
      prepack_pdb_title = target_ids[i] + "_tr_0001"

    # Generate a string of arguments from the case-specific variables
    s = Template( "python2 $rosetta../../tools/protein_tools/scripts/clean_pdb.py $title $chains")
    arguments = s.substitute( title=prepack_pdb_title, chains=partners, rosetta=config.rosetta_path )
    os.system( arguments )

def run_prepack_calc( energy_fxn, config, targets, test_name ): 
  """
    A function for generating prepacked structures for membrane protein-protein docking

    Arguments: 
      energy_fxn = energy function to use for calculations (typically, name of the weights file)
      config = path to benchmark, Rosetta executables
      targets = list of targets to use (homodimers, hetero-dimers unbound, hetero-dimers bound)
    """

  print("Generating prepacked structures for protein-protein complex set", targets)

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/structure/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
      temp = f.readlines()
      temp = [ x.strip() for x in temp ]
      target_ids = [ x.split("/")[0] for x in temp ]
      target_pdbs = [ x.split("/")[1] for x in temp ]

  # Generate path to executable
  executable = config.rosetta_path + "docking_prepack_protocol" + "." + config.platform + config.compiler + config.buildenv

  # Change directories to a data analysis dir
  outdir = config.benchmark_path + "/data/" + energy_fxn + "/protein-protein-docking/" 
  if ( not os.path.isdir( outdir ) ): 
    os.system( "mkdir " + outdir )
    os.chdir( outdir )

  outdir3 = outdir2 + "/" + targets
  if ( not os.path.isdir( outdir3 ) ): 
    os.system( "mkdir " + outdir3 )
    os.chdir( outdir3 )

  for i in range(0, len(target_ids)): 

    # Setup case-specific variables (pdbfile, spanfile)
    native = config.benchmark_path + "targets/structure/" + targets + "/" + target_ids[i] + "/" + target_pdbs[i]
    spanfile = config.benchmark_path + "targets/structure/" + targets + "/" + target_ids[i] + "/" + target_ids[i] + ".span"
    partners = target_pdbs[i].split("_")[1][0] + "_" + target_pdbs[i].split("_")[1][1]

    casedir = outdir3 + "/" + target_ids[i]
    if ( not os.path.isdir( casedir ) ): 
      os.system( "mkdir " + casedir )
      os.chdir( casedir )

    # Generate a string of arguments from the case-specific variables
    s = Template( " -in:file:s $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -docking:partners $partners -packing:pack_missing_sidechains 0 -nstruct 1 -out:path:all $outdir -mp:lipids:composition DLPC -mp:lipids:has_pore false" )
    arguments = s.substitute( native=native, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir, partners=partners )

    # Write jobfile and submit to the HPC
    print("Submitting docking calculations for case:", target_ids[i])
    jobname = target_ids[i] + "_prepack_dock"
    jobfile = hpc_util.make_jobfile( casedir, target_ids[i], executable, arguments )
    os.system( "bash " + jobfile )

def run_docking_calc( energy_fxn, config, targets, test_name, local_refine, transformed=False ): 
  """
  A function for generating prepacked structures for membrane protein-protein docking

  Arguments: 
    energy_fxn = energy function to use for calculations (typically, name of the weights file)
    config = path to benchmark, Rosetta executables
    targets = list of targets to use (homodimers, hetero-dimers unbound, hetero-dimers bound)
  """

  print("Generating prepacked structures for protein-protein complex set", targets)

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/structure/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
      temp = f.readlines()
      temp = [ x.strip() for x in temp ]
      target_ids = [ x.split("/")[0] for x in temp ]
      target_pdbs = [ x.split("/")[1] for x in temp ]

  # Generate path to executable
  executable = config.rosetta_path + "mpdocking" + "." + config.platform + config.compiler + config.buildenv

  # Change directories to a data analysis dir
  outdir = config.benchmark_path + "/data/" + energy_fxn + "/protein-protein-docking/" + targets + "/"
  if ( not os.path.isdir( outdir ) ): 
    sys.exit()

  for i in range(0, len(target_ids)): 

    # Setup case-specific variables (pdbfile, spanfile)
    native = config.benchmark_path + "targets/structure/" + targets + "/" + target_ids[i] + "/" + target_pdbs[i]
    spanfile = config.benchmark_path + "targets/structure/" + targets + "/" + target_ids[i] + "/" + target_ids[i] + ".span" 
    partners = target_pdbs[i].split("_")[1][0] + "_" + target_pdbs[i].split("_")[1][1]
    prepacked = config.benchmark_path + "data/protein-protein-docking/" + energy_fxn + "/" + targets + "/" + target_ids[i] + "/" + target_ids[i] + "_0001_" + target_pdbs[i].split("_")[1] + ".pdb"
    if ( transformed ): 
      prepacked = config.benchmark_path + "data/protein-protein-docking/" + energy_fxn + "/" + targets + "/" + target_ids[i] + "/" + target_ids[i] + "_tr_0001_" + target_pdbs[i].split("_")[1] + ".pdb"


    casedir = outdir + "/" + target_ids[i]
    if ( not os.path.isdir( casedir ) ): 
      sys.exit( "Prepack step was skipped! Exiting...")
    os.chdir( casedir )

    nmodels = 5000
    if ( local_refine ): 
      nmodels = 100

    # Generate a string of arguments from the case-specific variables
    s = Template( " -in:file:s $prepacked -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -run:multiple_processes_writing_to_one_directory -docking:partners $partners -docking:dock_pert 3 8 -packing:pack_missing_sidechains 0 -nstruct $nstruct -out:path:all $outdir -mp:lipids:composition DLPC -mp:lipids:has_pore false" )
    arguments = s.substitute( native=native, prepacked=prepacked, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir, partners=partners, nstruct=nmodels )

    # Add argument for local refine if applicable
    suffix = "_docking.sh"
    jobname = target_ids[i] + "_docking"
    if ( local_refine ): 
      arguments = arguments + " -docking_local_refine -out:prefix local_refine_"
      suffix = "_dock_local_refine.sh" 
      jobname = target_ids[i] + "_dock_local_refine"

    # Write jobfile and submit to the HPC
    print("Submitting docking calculations for case:", target_ids[i])
    jobfile = hpc_util.make_jobfile( casedir, target_ids[i], executable, arguments, suffix )
    hpc_util.submit_condor_job( casedir, jobname, jobfile, "", 100 )

def analyze_interfaces( energy_fxn, config, targets, test_name, local_refine, transformed=False ): 
  """
  A function for analyzing properties of the interfaces of docked models

  Arguments: 
    energy_fxn = energy function to use for calculations (typically, name of the weights file)
    config = path to benchmark, Rosetta executables
    targets = list of targets to use (homodimers, hetero-dimers unbound, hetero-dimers bound)
  """

  print("Analyzing interfaces protein-protein complex set", targets)

  # Read list of energy landscape test cases
  list_of_targets = config.benchmark_path + "targets/structure/" + targets + "/targets.list"
  with open( list_of_targets, 'rt' ) as f: 
      temp = f.readlines()
      temp = [ x.strip() for x in temp ]
      target_ids = [ x.split("/")[0] for x in temp ]
      target_pdbs = [ x.split("/")[1] for x in temp ]

  # Generate path to executable
  executable = config.rosetta_path + "rosetta_scripts" + "." + config.platform + config.compiler + config.buildenv
  xml_script = config.benchmark_path + "/src/xml/interface.xml"  

  # Change directories to a data analysis dir
  outdir = config.benchmark_path + "/data/" + energy_fxn + "/protein-protein-docking/"  + targets + "/"
  if ( not os.path.isdir( outdir ) ): 
    sys.exit()

  for i in range(0, len(target_ids)): 

    # Setup case-specific variables (pdbfile, spanfile)
    native = config.benchmark_path + "targets/structure/" + targets + "/" + target_ids[i] + "/" + target_pdbs[i]
    spanfile = config.benchmark_path + "targets/structure/" + targets + "/" + target_ids[i] + "/" + target_ids[i] + ".span" 
    partners = target_pdbs[i].split("_")[1][0] + "_" + target_pdbs[i].split("_")[1][1]

    casedir = outdir + "/" + target_ids[i]
    os.chdir( casedir )

    # Setup the scorefile name
    scorefile = casedir + "/" + target_ids[i] + "_interfaces.sc"
    if ( local_refine ): 
      scorefile = casedir + "/" + target_ids[i] + "interfaces_local.sc" 

    # Make list of models
    models_list = casedir + "/models.list" 
    structure_prefix =  casedir + "/" + target_ids[i] + "_0001_" + partners.split("_")[0] + partners.split("_")[1] + "*.pdb"
    if ( transformed ): 
      structure_prefix = casedir + "/" + target_ids[i] + "_tr_0001_" + target_pdbs[i].split("_")[1].split(".")[0] + "*.pdb"
    os.system( "ls " + structure_prefix + " > " + models_list )

    # Make a list of local refinement models
    local_models_list = casedir + "/models.list" 
    local_structure_prefix =  casedir + "/" + target_ids[i] + "_0001_" + partners.split("_")[0] + partners.split("_")[1] +  "*.pdb"
    if ( transformed ): 
      local_structure_prefix = casedir + "/" + target_ids[i] + "_tr_0001_" + target_pdbs[i].split("_")[1] + "*.pdb"
    os.system( "ls " + local_structure_prefix + " > " + local_models_list )
 
    # Generate a string of arguments from the case-specific variables
    s = Template( " -in:file:native $native -mp:setup:spanfiles $spanfile -score:weights $sfxn -out:file:scorefile $scorefile -out:file:score_only -parser:script_vars interface=$partners -packing:pack_missing_sidechains 0 -out:path:all $outdir -mp:lipids:composition DLPC -mp:lipids:has_pore false -parser:protocol $xml" )
    arguments = s.substitute( native=native, spanfile=spanfile, sfxn=energy_fxn, outdir=casedir, partners=partners, scorefile=scorefile, xml=xml_script )

    # Add argument for local refine if applicable
    suffix = "_interface_analysis.sh"
    jobname = target_ids[i] + "_interface_analysis"
    if ( local_refine ): 
      suffix = "_interface_analysis_local.sh" 
      jobname = target_ids[i] + "_interface_analysis_local"
      arguments = arguments + " -in:file:l " + local_models_list
    else: 
      arguments = arguments + " -in:file:l " + models_list

    # Write jobfile and submit to the HPC
    print("Submitting docking calculations for case:", target_ids[i])
    jobfile = hpc_util.make_jobfile( casedir, target_ids[i], executable, arguments, suffix )
    hpc_util.submit_condor_job( casedir, jobname, jobfile, "", 1 )
