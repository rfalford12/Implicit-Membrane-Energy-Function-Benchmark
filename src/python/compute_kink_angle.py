"""
MeasureDegreeKink:
1.   Method to find the dotproduct of vectors
2.   Method to find the length of a vector
3.   Method to find the angle between two vectors
4.   Method to find the kink point or hinge residue within a helix
5.   Method to find the backbone coordinates of the helix
6.   Method to find the angle of the kink
"""
import rosetta
from pyrosetta import*
init()
from pyrosetta import toolbox
from pyrosetta.rosetta.protocols.membrane import*
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import numpy as np
import math
class MeasureDegreeKink():
    """Measures the degree to which a residue kinks"""
    def __init__(self, file):
        """Construct the class"""
        self.file=file
        self.pose=pose_from_file(self.file)
    def dotproduct(self,v1,v2):
        """ Dot product of two vectors"""
        # 1. Input vectors to find the dot product
        return sum((a*b) for a, b in zip(v1, v2))
    def length(self,v):
        """ Length of a vector """
        # 2. Input a vector and return its length
        return math.sqrt(MeasureDegreeKink(self.file).dotproduct(v, v))
    def angle(self,v1, v2):
        """ Angle between two vectors """
        # 3. Take two vectors and output the angle between them
        return math.acos(MeasureDegreeKink(self.file).dotproduct(v1, v2) / (MeasureDegreeKink(self.file).length(v1) * MeasureDegreeKink(self.file).length(v2)))
    def kink_point(self):
        """ Determines the kink point within a helix """
        # 4. Determine the hinge residue for kinking within a helix
        #    by finding which residue has the greatest phi and psi
        #    change from the ideal dihedral angles of an alpha helix
        phi_ave = -64.0
        psi_ave = -41.0
        set=0
        res_num=1
        end=self.pose.total_residue()+1
        for x in xrange(1,end):
            phi=self.pose.phi(x)
            psi=self.pose.psi(x)
            dif_phi = math.fabs(phi_ave-phi)
            dif_psi = math.fabs(psi_ave-psi)
            tot=dif_phi+dif_psi
            #print tot
            if tot > set and x >= 5 and x <= self.pose.total_residue()-4:
                res_num = x
                set = dif_phi+dif_psi
            else:
                set=set
                res_num=res_num
        return res_num
        #print res_num
    def direction(self):
        """ Find the direction of PCA vectors """
        # 4. Using vectors from the start to kink coords or the kink to end coords
        #    the angle can be found to determine if the PCA vectors 1 and 2 are
        #    pointing in the right direction
        coord = MeasureDegreeKink(self.file).count(0,self.pose.total_residue())
        kink_point = MeasureDegreeKink(self.file).kink_point()
        point1 = np.array(coord[0])
        point2 = np.array(coord[kink_point-1])
        point3 = np.array(coord[-1])
        vector12 = point2-point1
        vector23 = point3-point2
        coord1 = coord[0:kink_point+1]
        coord2 = coord[kink_point:len(coord)+1]
        pca1 = PCA(n_components=1)
        pca2 = PCA(n_components=1)
        pca1.fit(coord1)
        pca2.fit(coord2)
        Vpca1 = pca1.components_[0]
        Vpca2 = pca2.components_[0]
        diff1 = MeasureDegreeKink(self.fragfile).angle(vector12,Vpca1)*180/3.14
        diff2 = MeasureDegreeKink(self.fragfile).angle(vector23,Vpca2)*180/3.14
        diffs=[diff1,diff2]
        return diffs
    def count(self, start, end):
        """ Find the coordinates for the backbone atoms of the helix"""
        # 5. Find the backbone coordinates of the helix by looping through
        #    each residue and finding the coordinates of the Calpha, the N
        #    and C of each residue
        coord=[]
        for x in xrange(start,end):
            CA=self.pose.residue(x).xyz("CA")
            N=self.pose.residue(x).xyz("N")
            C=self.pose.residue(x).xyz("C")
            x1=[CA[0],CA[1],CA[2]]
            x2=[N[0],N[1],N[2]]
            x3=[C[0],C[1],C[2]]
            coord.append(x1)
            coord.append(x2)
            coord.append(x3)
        coord=np.array(coord)
        return coord
    def kink_angle(self):
        """ Determines the kink angle of a helix by using PCA"""
        # 6. Use Principle Component Analysis (PCA)to find a vector
        #    that points in the direction of the helix before the
        #    hinge point, and a vector that points in the direction
        #    of the helix after the hinge point. Using the two vectors
        #    find the angle between the two vectors
        #print "this"
        start = MeasureDegreeKink(self.file).kink_point()
        #print start
        end = MeasureDegreeKink(self.file).kink_point()+1
        #print end
        coord1 = MeasureDegreeKink(self.file).count(1,end)
        #print coord1
        coord2 = MeasureDegreeKink(self.file).count(start,self.pose.total_residue())
        pca1 = PCA(n_components=1)
        pca2 = PCA(n_components=1)
        pca1.fit(coord1)
        pca2.fit(coord2)
        vector1 = pca1.components_[0]
        vector2 = pca2.components_[0]
        #print vector1
        #print vector2
        angle=MeasureDegreeKink(self.file).angle(vector1,vector2)*180/3.14
        #print angle
        if angle >= 90:
            angle=math.fabs(180-angle)
        else:
            angle=angle
        print "Angle (degree): "+ str(angle)
MeasureDegreeKink("/Users/Maestro/apps/PyRosetta4/3f5wEndHelix.pdb").kink_angle()
