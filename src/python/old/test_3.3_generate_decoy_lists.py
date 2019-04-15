#!/usr/bin/env python
# @file: test_3.3_generate_decoy_lists.py
# @brief: Take a master list of decoys and subdivide into 50 smaller lists for large decoy discrimination jobs
# @author: Rebecca F. Alford (ralford3@jhu.edu)

import sys, os
from string import Template
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def chunks(l, n):
    n = max(1, n)
    return (l[i:i+n] for i in xrange(0, len(l), n))

def main( args ): 

    parser = OptionParser(usage="usage %prog --decoy_list AcrB_models.list --path /home/ralford3" )
    parser.set_description(main.__doc__)

    parser.add_option('--decoy_list', '-l', 
        action="store",
        help="List of decoys to subdivide", )

    parser.add_option( "--path", '-p', 
        action="store", 
        help="Path to prepend each model", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that the required options have been provided
    if ( not Options.decoy_list or not Options.path ): 
        print "Missing required options --path or --decoy_list! Exiting..."
        sys.exit()

    # read the file into a list wthout the old paths
    with open( Options.decoy_list, 'rb' ) as f: 
        decoys = f.readlines()
        decoys = [ x.strip() for x in decoys ]
        decoys = [ x.split("/") for x in decoys ]
        decoys = [ x[-1] for x in decoys ]

    # divide the list into 50 chunks
    subdivide = list(chunks( decoys, 50 ))
    
    # Write the output files
    for i in xrange(1, 51): 
        outfile = "models." + str(i) + ".list" 
        with open( outfile, 'a' ) as f: 
            for item in subdivide[i]: 
                f.write( Options.path + "/" + item + "\n" )
            f.close()

if __name__ == "__main__" : main(sys.argv)


