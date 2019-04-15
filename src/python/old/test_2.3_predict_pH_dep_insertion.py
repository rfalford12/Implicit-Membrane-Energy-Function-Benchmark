# @file: predict_ddG_of_insertion.py
# @brief: Calculate the ddG of inserting a peptide into the bilayer from solution
# @author: Rebecca Alford (ralford3@jhu.edu)

import sys, os
import commands
import random
from string import Template

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main( args ):

    # Read options from the commandline
    parser = OptionParser( usage="usage: %prog --rscript landscape.xml --pdb test.pdb --span span.pdb --energy_fxn mpframework_fa_2007" )
    parser.set_description(main.__doc__)

    parser.add_option('--infile_pH4', '-a',
        action="store", 
        help="Input energy landscape at pH 4",)

    parser.add_option('--infile_pH7', '-b', 
        action="store", 
        help="Input energy landscape at pH 7",)

    parser.add_option('--outfile', '-o',
        action="store",
        help="Output file containing ddg data",)

    parser.add_option('--case', '-c', 
        action="store", 
        help="Name of test case",)

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    case = Options.case

    # Generate a 2D energy landscpae at both pH = 4 and pH = 7
    pH4_zcoords, pH4_angles, pH4_scores = generate_pH_dependent_landscape( Options.infile_pH4 )
    pH7_zcoords, pH7_angles, pH7_scores = generate_pH_dependent_landscape( Options.infile_pH7 )

    # Get the score of the non-inserted helix
    max_index = pH4_zcoords.index( max(pH4_zcoords) )
    non_inserted_score = pH4_scores[ max_index ]

    # Get the lowest energy coordinate of the pH7 helix
    min_index = pH7_scores.index( min(pH7_scores) )
    lowest_zcoord = pH7_zcoords[ min_index ]
    lowest_angle = pH7_angles[ min_index ]
    lowest_score = pH7_scores[ min_index ]
    ddG_from_lowest = lowest_score - non_inserted_score

    # Get the score of the membrane centered pose ("fully inserted")
    rounded_zcoord = [ round(x) for x in pH7_zcoords ]
    inserted_index = rounded_zcoord.index( 0 )
    inserted_zcoord = pH7_zcoords[ inserted_index ]
    inserted_angle = pH7_angles[ inserted_index ]
    inserted_score = pH7_scores[ inserted_index ]
    ddG_from_inserted = inserted_score - non_inserted_score

    with open( Options.outfile, 'a' ) as f:
        s = Template( "$case $non_inserted_score $lowest_zcoord $lowest_angle $lowest_score $ddG_from_lowest $inserted_zcoord $inserted_angle $inserted_score $ddG_from_inserted\n" )
        outstr = s.substitute( case=case, non_inserted_score=non_inserted_score, lowest_zcoord=lowest_zcoord, lowest_angle=lowest_angle, lowest_score=lowest_score, ddG_from_lowest=ddG_from_lowest, inserted_zcoord=inserted_zcoord, inserted_angle=inserted_angle, inserted_score=inserted_score, ddG_from_inserted=ddG_from_inserted )
        f.write( outstr )

def generate_pH_dependent_landscape( pH_output ):

    # Read file into 3D data structure
    with open( pH_output, 'r' ) as f:
        content = f.readlines()
    content = [ x.strip() for x in content ]
    content = [ x.split(' ') for x in content ]

    # Read lines into an array of data triplets
    zcoords = []
    angles = []
    scores = []
    for x in content:
        if ( x[0] != "zcoord" ):
            zcoords.append( float(x[0]) )
            angles.append( float(x[1]) )
            scores.append( float(x[2]) )

    return zcoords, angles, scores

if __name__ == "__main__" : main(sys.argv)
