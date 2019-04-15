#!/usr/bin/env python
#@file: test_3.3_process_decoy_disc_data.py
#@brief: Process decoy discrimination data for large and small sets
#@author: Rebecca F. Alford (ralford3@jhu.edu)

import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

### Global Data ###
boincdir = "/home/ralford/apps/Rosetta/bakerlab_scripts/boinc/"
benchmark = "/home/ralford/membrane_efxn_benchmark/"
rosettadir = "/home/ralford/apps/Rosetta/main/source/bin/"
buildenv = "release"
compiler = "gcc"

def main( args ): 

    parser = OptionParser(usage="usage %prog --energy_fxn membrane_01" )
    parser.set_description(main.__doc__)

    parser.add_option('--energy_fxn', '-e',
        action="store",
        help="Name of energy function dataset to post process", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    if ( not Options.energy_fxn ):
        print "Missing required option --energy_fxn! Exiting..."
        sys.exit()

    process_decoy_disc_results( benchmark, Options.energy_fxn )

def process_decoy_disc_results( benchmark, energy_fxn ):

    print "Processing decoy discrimination results generated with the energy function", energy_fxn

    # Change directories into a data analysis dir
    path_to_test = benchmark + "inputs/test_3.3_decoy_discrimination/"
    outdir = benchmark + "data/" + energy_fxn + "/test_3.3_decoy_discrimination"
    os.chdir( outdir )

    # Setup global rescoring variables
    executable = rosettadir + "score_jd2.linux" + compiler + buildenv
    r = Template( " -overwrite -in:file:l $refined_list -in:file:native $native -in:membrane -mp:setup:spanfiles $spanfile -out:file:scorefile $output_sc -score:weights $energy_fxn -mp:lipids:has_pore false" )

    # Process Small-Decoy Dataset (From Dutagaci et al.)
    list_of_test_cases = path_to_test  + "/dutagaci-set/decoy_sets.list"
    with open( list_of_test_cases, 'rb' ) as f: 
        test_cases = f.readlines()
    test_cases = [ x.strip() for x in test_cases ]

    for case in test_cases: 
        casedir = outdir + "/dutagaci-set/" + case 
        os.chdir( casedir )

        #Setup variables for rescoring
        os.system( "ls -d " + casedir + "/*_0001.pdb > models_list" )
        models_list = casedir + "/models_list"
        native = path_to_test + "/dutagaci-set/" + case + "/" + case + "_native.pdb"
        spanfile = path_to_test + "/dutagaci-set/" + case + "/" + case + ".span"
        output_sc = casedir + "/" + case + "_rescored.sc"
        arguments0 = r.substitute( refined_list=models_list, native=native, spanfile=spanfile, output_sc=output_sc, energy_fxn=energy_fxn )
        
        # Rescore all of the refined decoys to calculate an rms relative to the native
        rescore_cmd = executable + " " + arguments0
        os.system( rescore_cmd )


        # Run decoy discrimination on each parsed scorefile (here, its by Irms and I_sc)
        scoring_script = boincdir + "score_energy_landscape.py"
        output_discfile = casedir + "/" + case + ".disc"
        t = Template( " -terms rms total_score -abinitio_scorefile $docking_sc > $docking_disc" )
        arguments2 = t.substitute( docking_sc=output_sc, docking_disc=output_discfile )
        os.system( "python " + scoring_script + arguments2 )

    ## Process Large-Decoy Dataset (From Yarov-Yaraovy et al.)
   # list_of_test_cases = path_to_test  + "/yarov-yaravoy-set/decoy_sets.list"
   # with open( list_of_test_cases, 'rb' ) as f:
   #     test_cases = f.readlines()
   # test_cases = [ x.strip() for x in test_cases ]

#    for case in test_cases:
#        casedir = outdir + "/yarov-yaravoy-set/" + case
  #      os.chdir( casedir )

        # Setup variables for rescoring                                                                                              
#        os.system( "ls -d " + casedir + "/*_0001.pdb > models_list" )
#        models_list = casedir + "/models_list"
#        native = path_to_test + "/yarov-yaravoy-set/" + case + "/" + case + "_native.pdb"
#        spanfile = path_to_test + "/yarov-yaravoy-set/" + case + "/" + case + ".span"
#        output_sc = casedir + "/" + case + "_rescored.sc"
#        arguments0 = r.substitute( refined_list=models_list, native=native, spanfile=spanfile, output_sc=output_sc, energy_fxn=energy_fxn )
        
        # Rescore all of the refined decoys to calculate an rms relative to the native                         
 #       rescore_cmd = executable + " " + arguments0
 #       os.system( rescore_cmd )

        # Run decoy discrimination on each parsed scorefile (here, its by Irms and I_sc)
  #      scoring_script = boincdir + "score_energy_landscape.py"
  #      output_discfile = casedir + "/" + case + ".disc"
  #      t = Template( " -terms rms total_score -abinitio_scorefile $docking_sc > $docking_disc" )
  #      arguments2 = t.substitute( docking_sc=output_sc, docking_disc=output_discfile )
  #      os.system( "python " + scoring_script + arguments2 )

        
if __name__ == "__main__" : main(sys.argv)
