#@file: parse_scorefile.py
#@brief: Given a docking scorefile, only print the header and lines containing 40 elements
#author: Rebecca F. Alford (ralford3@jhu.edu)

import sys, os
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

from pyrosetta import *
from pyrosetta.teaching import *

## @brief: Parse scorefile by number of expected elements
def main( args ):

    parser = OptionParser( usage="usage: %prog --models docking_list --output 1afo.labels.txt" )
    parser.set_description(main.__doc__)

    #input options
    parser.add_option('--models', '-m',
        action="store",
        help="Input list of docking models",)

    parser.add_option('--spanfile', '-s',
        action="store",
        help="Input list of "
    )

    parser.add_option('--output', '-o',
        action="store",
        help="Name of output file contianing per-model insertion labels",)

    #parse options
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    #error checking
    if ( not Options.models or not Options.spanfile or not Options.output ):
        print "Missing required option --models, --spanfile, or --output! Exiting..."
        sys.exit()

    if ( not os.path.isfile( Options.models ) or not os.path.isfile( Options.spanfile ) ):
        print "Filepaths for models list or spanfile do not exist! Exiting..."
        sys.exit()

    init( extra_options="-out:mute all -in:ignore_unrecognized_res")

    # read list of models line by line
    with open( Options.models, 'rb' ) as f:
        content = f.readlines()
    models = [ x.strip() for x in content ]

    # Open an output file for Writing
    f = open( Options.output, 'a' )
    f.write( "Model Inserted\n" )

    for pdb in models:

        # Read in docked model from file
        pose = pose_from_pdb( pdb )
        add_memb = rosetta.protocols.membrane.AddMembraneMover( Options.spanfile )
        add_memb.apply( pose )

        # make a vector1 of residue z coordinates
        res_z_coord = rosetta.utility.vector1_double()
        for i in range(1, pose.total_residue()+1):
            res_z_coord.append( pose.residue(i).atom("CA").xyz().z )

        # Determine if the pose is inserted into the membrane
        in_membrane = True
        spanning_topology = pose.conformation().membrane_info().spanning_topology()
        for i in range(1, spanning_topology.nspans() ):
            if ( not spanning_topology.spanning( res_z_coord, spanning_topology.get_spans()[i] ) ):
                in_membrane = False

        print pdb.split("/")[-1], in_membrane

        # Write the result to a file
        f.write( pdb.split("/")[-1] + " " + str(in_membrane) )

    f.close()


if __name__ == "__main__" : main(sys.argv)
