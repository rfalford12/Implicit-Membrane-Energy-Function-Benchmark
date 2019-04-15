#!/usr/bin/env python
#@file: test_3.2_process_docking_data.py
#@brief: Parse the docking scorefile and calculate decoy discrimination
#@author: Rebecca F. Alford (ralford3@jhu.edu)

import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

### Global Data ###
boincdir = "/Volumes/ralford/apps/Rosetta/bakerlab_scripts/boinc/"
benchmark = "/Volumes/ralford/membrane_efxn_benchmark/"

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

    process_docking_results( benchmark, Options.energy_fxn )


def process_docking_results( benchmark, energy_fxn ):

    print "Processing docking results generated with the energy function", energy_fxn

    # Change directories into a data analysis dir
    outdir = benchmark + "data/" + energy_fxn + "/test_3.2_docking/"
    os.chdir( outdir )

    # Process the small homodimer dataset
    small_homodimers_list = benchmark + "inputs/test_3.2_docking/small-homodimer-set/dimers.list"
    with open( small_homodimers_list, 'rb' ) as f: 
        small_homodimers = f.readlines()
        small_homodimers = [ x.strip() for x in small_homodimers ]
    
    for case in small_homodimers: 

        print "Estimating decoy discrimination for small homodimer case", case

        # Make a path for the homodimer case
        casedir = outdir + "small-homodimer-set/" + case 
        os.chdir( casedir )

        # Parse each scorefile to exclude incomplete lines
        input_path = casedir + "/score.sc"
        if ( os.path.isfile( input_path ) ): 
            output_path = casedir + "/" + case + "_score.sc"
            exclude_incompletely_scored_decoys( input_path, output_path )

            # Run decoy discrimination on each parsed scorefile (here, its by Irms and I_sc)
            scoring_script = boincdir + "score_energy_landscape.py"
            output_discfile = casedir + "/" + case + ".tdisc"
            t = Template( " -terms rms total_score -abinitio_scorefile $docking_sc > $docking_disc" )
            arguments2 = t.substitute( docking_sc=output_path, docking_disc=output_discfile )
            os.system( "python " + scoring_script + arguments2 )
        else: 
            "Skipping this case: ", case

    # Now repeat for the large homodimer dataset
    large_homodimers_list = benchmark + "inputs/test_3.2_docking/large-homodimer-set/dimers.list"
    with open( large_homodimers_list, 'rb' ) as f: 
        large_homodimers = f.readlines()
        large_homodimers = [ x.strip() for x in large_homodimers ]

    for case in large_homodimers: 

        print "Estimating decoy discrimination for large homodimer case", case

        # Make a path for the homodimer case
        casedir = outdir + "large-homodimer-set/" + case 
        os.chdir( casedir )

        # Parse each scorefile to exclude incomplete lines
        input_path = casedir + "/score.sc"
        if ( os.path.isfile( input_path ) ): 
            output_path = casedir + "/" + case + "_score.sc"
            exclude_incompletely_scored_decoys( input_path, output_path )

            # Run decoy discrimination on each parsed scorefile (here, its by Irms and I_sc)
            scoring_script = boincdir + "score_energy_landscape.py"
            output_discfile = casedir + "/" + case + ".tdisc"
            t = Template( " -terms rms total_score -abinitio_scorefile $docking_sc > $docking_disc" )
            arguments2 = t.substitute( docking_sc=output_path, docking_disc=output_discfile )
            os.system( "python " + scoring_script + arguments2 )
        else: 
            "Skipping this case: ", case

def exclude_incompletely_scored_decoys( infile, outfile ): 

    # read the file line by line
    with open( infile, 'rb' ) as f: 
        content = f.readlines()
        content = [ x.strip() for x in content ]

    # determine the expected number of elements
    n_elements = len(content[1].split())

    # only accept elements with the desired number of lines
    with open( outfile, 'w' ) as f: 
        for line in content: 
            if ( len(line.split()) == n_elements ): 
                f.write( line + "\n" )

if __name__ == "__main__" : main(sys.argv)
