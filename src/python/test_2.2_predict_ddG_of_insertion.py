# @file: predict_ddG_of_insertion.py
# @brief: Calculate the ddG of inserting a peptide into the bilayer from solution
# @author: Rebecca Alford (ralford3@jhu.edu)

#from pyrosetta import *
from string import Template

import sys, os
import commands
import random

from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def main( args ):

    # Read options from the commandline
    parser = OptionParser( usage="usage: %prog" )
    parser.set_description(main.__doc__)

    parser.add_option('--infile', '-i', 
        action="store", 
        help="Input energy landscape file",)

    parser.add_option('--case', '-c', 
        action="store", 
        help="Name of test case",)

    parser.add_option('--outfile', '-o',
        action="store",
        help="Output file containing ddg data",)

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Read file into 3D data structure
    landscapefile = Options.infile
    case = Options.case
    with open( landscapefile, 'r' ) as f:
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

    # Get the energy of the non-inserted pose
    max_index = zcoords.index( max(zcoords) )
    non_inserted_score = scores[ max_index ]

    # Get the lowest energy coordinate
    min_index = scores.index( min(scores) )
    lowest_zcoord = zcoords[ min_index ]
    lowest_angle = angles[ min_index ]
    lowest_score = scores[ min_index ]
    ddG_from_lowest = lowest_score - non_inserted_score

    # Get the score of the membrane centered pose ("fully inserted")
    rounded_zcoord = [ round(x) for x in zcoords ]
    inserted_index = rounded_zcoord.index( 0 )
    inserted_zcoord = zcoords[ inserted_index ]
    inserted_angle = angles[ inserted_index ]
    inserted_score = scores[ inserted_index ]
    ddG_from_inserted = inserted_score - non_inserted_score

    with open( Options.outfile, 'a' ) as f:
        s = Template( "$case $non_inserted_score $lowest_zcoord $lowest_angle $lowest_score $ddG_from_lowest $inserted_zcoord $inserted_angle $inserted_score $ddG_from_inserted\n" )
        outstr = s.substitute( case=case, non_inserted_score=non_inserted_score, lowest_zcoord=lowest_zcoord, lowest_angle=lowest_angle, lowest_score=lowest_score, ddG_from_lowest=ddG_from_lowest, inserted_zcoord=inserted_zcoord, inserted_angle=inserted_angle, inserted_score=inserted_score, ddG_from_inserted=ddG_from_inserted )
        f.write( outstr )

if __name__ == "__main__" : main(sys.argv)
