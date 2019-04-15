#@file: parse_scorefile.py
#@brief: Given a docking scorefile, only print the header and lines containing 40 elements
#author: Rebecca F. Alford (ralford3@jhu.edu)

import sys, os
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

## @brief: Parse scorefile by number of expected elements
def main( args ):

    parser = OptionParser( usage="usage: %prog --input score.sc --num_elements 40 ")
    parser.set_description(main.__doc__)

    #input options
    parser.add_option('--input', '-i',
        action="store",
        help="Input scorefile",)

    parser.add_option('--num_elements', '-n',
        action="store",
        help="Number of elements expected per line",
    )
    parser.add_option('--output', '-o',
        action="store",
        help="Name of output scorefile"
    )

    #parse options
    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    #error checking
    if ( not Options.input or not Options.num_elements ):
        print "Missing required arguments --input or --num_elements! Exiting..."
        sys.exit()

    if ( not os.path.isfile( Options.input ) ):
        print "Incorrect path to user specified input file"
        sys.exit()

    num_elements = int( Options.num_elements )
    with open( Options.input, 'r' ) as f:
        content = f.readlines()
    content = [ x.strip() for x in content ]

    with open( Options.output, 'w' ) as f:
        for line in content:
            if ( len( line.split() ) == int(Options.num_elements) ):
                f.write( line + "\n" )

if __name__ == "__main__" : main(sys.argv)
