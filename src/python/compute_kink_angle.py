#!/usr/bin/env python
# @file: compute_kink_angle.py
# @author: Rebecca F. Alford (ralford3@jhu.edu), adapted from B. Lasher
# #brief: Compute the degree of kinking for a particular helix

from pyrosetta import*
from pyrosetta import toolbox
from pyrosetta.rosetta.protocols.membrane import*
import rosetta

from sklearn.decomposition import PCA
from optparse import OptionParser, IndentedHelpFormatter
import numpy as np
import math

init()

def get_dot_product(v1,v2):
    """
    Compute the dot product of two vectors
      - inputs: two vectors, v1 and v2
    """
    return sum((a*b) for a, b in zip(v1, v2))

def get_length(v):
    """
    Compute the length of a vector
      - inputs: one vector (v)
    """
    return math.sqrt(get_dot_product(v, v))

def get_angle(v1, v2):
    """
    Compute the angle between two vectors
      - inputs: two vectors, v1 and v2
    """
    return math.acos(get_dot_product(v1, v2) / (get_length(v1) * get_length(v2)))

def get_kink_point(pose, helix_start, helix_end):
    """
    Compute the hinge point within a helix 
    @details Determine the hing residue for a helix by finding
        the residue with the largest phi and psi change from 
        ideal dihedral angles
    """

    # Ideal phi and psi angles 
    phi_ave = -64.0
    psi_ave = -41.0

    ## not sure what this does
    current_diff = 0
    hinge_residue = 1

    for x in xrange(helix_start, helix_end+1): 

        # Calculate phi and psi of the current residue
        curr_phi = pose.phi(x)
        curr_psi = pose.psi(x)

        # Compute teh difference between the ideal and current phi/psi
        diff_phi = math.fabs(phi_ave-curr_phi)
        diff_psi = math.fabs(psi_ave-curr_psi)
        total_diff = diff_phi + diff_psi 

        # todo - adjust this line
        if (total_diff > current_diff) and (x >= 5) and (x <=pose.total_residue()-4): 
            hinge_residue = x
            current_diff = total_diff

    return hinge_residue

def get_coords(pose, start, end):
    """
    Find the coordinates for the backbone atom for a specific residue
    @details Get the backbone coordinates by looping through each residue
        and storing the Calpha, N, and C coordinates of each residue
    """

    # List to store coordiantes
    coord_list = []

    # Loop through the residues and store the coordinates
    for x in xrange(start, end):

        # Grab the xyz coordiantes from the pose
        CA_xyz = pose.residue(x).xyz("CA")
        N_xyz = pose.residue(x).xyz("N")
        C_xyz = pose.residue(x).xyz("C")

        # Store the coordinates in individual vectors
        x1 = [ CA_xyz[0], CA_xyz[1], CA_xyz[2] ]
        x2 = [ N_xyz[0], N_xyz[1], N_xyz[2] ]
        x3 = [ C_xyz[0], C_xyz[1], C_xyz[2] ]

        # Append coordinates to the list 
        coord_list.append(x1)
        coord_list.append(x2)
        coord_list.append(x3)

    coord_list = np.array(coord_list)
    return coord_list

def direction(pose, helix_start, helix_end): 
    """
    Find the direction of PCA vectors
    @details Using vectors that represent the hinge to start or end residues, 
        the angle can be found to determine of the PCA vectors 1 and 2 are
        pointing in the right direction
    """

    ### this is definitely the messiest

    # Get list of helix coordinates
    coord_list = get_coords(pose, helix_start, helix_end)
    kink_point = get_kink_point(pose, helix_start, helix_end)

    # Calculate the hinge point and the +1 and -1 positions
    point1 = np.array(coord_list[0])
    point2 = np.array(coord_list[kink_point-1])
    point3 = np.array(coord_list[-1])
    vector12 = point2-point1
    vector23 = point3-point2
    coord1 = coord_list[0:kink_point+1]
    coord2 = coord_list[kink_point:len(coord)+1]
    pca1 = PCA(n_components=1)
    pca2 = PCA(n_components=1)
    pca1.fit(coord1)
    pca2.fit(coord2)
    Vpca1 = pca1.components_[0]
    Vpca2 = pca2.components_[0]
    diff1 = angle(vector12,Vpca1)*180/3.14
    diff2 = angle(vector23,Vpca2)*180/3.14
    diffs=[diff1,diff2]
    return diffs

def compute_kink_angle(pose, helix_no): 
    """
    Calculate the helix kink angle using PCA
    @details Using principle component analysis (PCA), find a vector
        that points in the direction of the helix before and after 
        the hinge point. Then, find the angle between the two PCA vectors. 
    """

    # Calculate the start and end position of the helix
    topology = pose.conformation().membrane_info().spanning_topology()
    span = topology.span( helix_no )
    helix_start = span.start()
    helix_end = span.end()

    # Compute the kink point start/end coordinates
    hinge_start = get_kink_point(pose, helix_start, helix_end)
    hinge_end = helix_start+1

    # Error checking
    if ( helix_end == hinge_start ):
        if ( helix_end != pose.total_residue() ):  
            helix_end = helix_end + 1
        else: 
            return 0
    if ( helix_start == hinge_end ):
        if ( helix_start != 1 ):
            helix_start = helix_start - 1 
        else: 
            return 0

    coord1 = get_coords(pose, helix_start, hinge_end)
    coord2 = get_coords(pose, hinge_start, helix_end)

    # Setup and compute principal components
    pca1 = PCA(n_components=1)
    pca2 = PCA(n_components=1)
    pca1.fit(coord1)
    pca2.fit(coord2)
    vector1 = pca1.components_[0]
    vector2 = pca2.components_[0]

    # Compute the angle between the two vectors
    angle = get_angle(vector1,vector2)*180/3.14
    if angle >= 90:
        angle = math.fabs(180-angle)
    else:
        angle = angle
    
    return(angle)

def main( args ): 

    # Read options from the command line
    parser = OptionParser(usage="usage %prog --pdblist sample.list --spanfile sample.span" )
    parser.add_option( '--pdblist', '-l', action="store", help="List of PDB files" )
    parser.add_option( '--spanfile', '-s', action="store", help="Corresponding span file" )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Read list of pdbfiles
    with open( Options.pdblist, 'rt' ) as f: 
        contents = f.readlines()
        pdbfiles = [ x.strip() for x in contents ]

    # Open output file
    with open( "output_kinks_native.txt", 'wt' ) as g: 
        for pdb in pdbfiles: 
            # Read pdb from file
            if ( not os.path.isfile( pdb ) ): 
                print("PDB file does not exist!")
                sys.exit()
            print("Calculating kink angle for model", pdb)
            pose = pose_from_file( pdb ) 
            g.write( pdb + " " )

            # Add membrane to pdb file
            add_memb = rosetta.protocols.membrane.AddMembraneMover( Options.spanfile )
            add_memb.apply( pose )

            # For each span, calculate the kinking angles
            topology = pose.conformation().membrane_info().spanning_topology()
            nspans = topology.nspans()
            for i in range(1, nspans+1):
                kink_angle = compute_kink_angle(pose, i)
                kink_angle = round(kink_angle, 3)
                if ( i != nspans ): 
                    g.write( str(kink_angle) + " " )
                else: 
                    g.write( str(kink_angle) + "\n" )

if __name__ == "__main__" : main(sys.argv)
