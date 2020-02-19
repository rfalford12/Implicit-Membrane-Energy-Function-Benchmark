import sys
import os
import math
import subprocess as sub
from glob import glob
##for writing a csv file
import csv
import numpy as np
import argparse

from prosci.util.pdb import Pdb #parses pdbfile
from prosci.util.residue import ResidueList
from prosci.util.ali import Ali #parses tem_file
from prosci.util.pdb3d import fit_line


#############################################################
# Some definitions
#############################################################
def isGap(c):
  return c in "-/"

def isLoop(c):
  return c in "CPE"
  #return c != "H"

def angle(crds1, crds2):
  return math.acos(np.dot(crds1, crds2)/(np.linalg.norm(crds1)*
                                            np.linalg.norm(crds2)))

def separation_distance(crds1,crds2):
  return np.linalg.norm(crds1-crds2)
  
def dist(a, b):
    return math.sqrt((a.x-b.x)**2 + (a.y-b.y)**2 + (a.z-b.z)**2)

def wobble2(C, R_nfit, n_unitvector, c_unitvector):
  #point C is the last C alpha atom in the n_fragment
  #calculates the 'wobble' angle of a helix kink
  #180 is outside the kink!
  
  #point X where the plane, perpendicular to n_unitvector, which contains
  # C, intersects with n_vector is given by n_linepoints[1] + n_unit_vector * 
  #distance(C, n_lineopoints[1]) * sin (alpha-pi/2)
  point_X=R_nfit+n_unitvector*(np.dot(n_unitvector,(C-R_nfit)))
  
  #then the direction of the projection of the c unit vector onto the plane is 
  #given by
  XP=np.cross(n_unitvector, np.cross(n_unitvector,c_unitvector))
  #and the vector between the c alpha atom, and point X is:
  XC=(point_X - C)
  
  #then the 2nd angle is given by
  wobble_ang=math.degrees(angle(XP,XC))
  
  #but this is only in the range 0 to pi, and in fact, we want the 0 to 2 pi to
  #fully describe this angle
  #the angle between XP x XC and n_unitvector will be either 0 or 180, which 
  #determines whether we are in 0 - pi or pi - 2pi 
  if np.dot(n_unitvector, np.cross(XP,XC))>0:
    wobble_ang = 360-wobble_ang  
  
  return wobble_ang

def chop_helices(layers): # now split by how far from the H layer
  #take the first and last instances of not N
  for i in xrange(len(layers)):
    if layers[i] != 'N':
      membrane_start = i
      break
  for i in xrange(len(layers)-1, -1,-1): #and the last instance of not N
    if layers[i] != 'N':
      membrane_end = i
      break      
  helix_start = max(0,membrane_start-5) #extend 5 either way
  
  helix_end = min(len(layers),membrane_end+5)    
  
  return [helix_start,helix_end, membrane_start, -1*membrane_end] #done

def find_helices(pdb,sstruc, sequence, layers, soluble, max_loop_length):
  sstruc_bounds = []
  for i in xrange(1,len(sstruc)):
    if not isGap(sequence[i]):
      if not isLoop(sstruc[i]):
        # Either: There was a break in the protein chain
        if dist(pdb[i-1].CA, pdb[i].CA) > 4:
          sstruc_bounds.append([i, i+1])
        
        #elif  abs(pdb[i-1].CA.ires - pdb[i].CA.ires) != 1: #not sequential residue numbers:
          #sstruc_bounds.append([i, i+1]) 
        
        # Or: We're starting a new secondary structure element
        # (not part for first iteration)
        elif not sstruc_bounds or sstruc_bounds[-1][1] < i:
          
          # Jumping over a small loop 
          #again, first part is for when no sstruc exists
          if sstruc_bounds and \
          i-sstruc_bounds[-1][1] < max_loop_length:
            # We've just jumped a tiny loop. Merge this SSE with the 
            #previous one and delete the tiny loop.
            sstruc_bounds[-1][1] = i+1
            
          # Really a new secondary structure element
          else:
            sstruc_bounds.append([i, i+1])
        
        # Or: Extending SSE
        else:
          sstruc_bounds[-1][1] = i+1
  # Find TM helices - pick ones only with part in the tail region in the membrane
  tm_helices = []
  if not soluble: 
    for start, end in sstruc_bounds:
      is_tm_helix = "T" in layers[start:end] and "H" in sstruc[start:end]
      if is_tm_helix:
        new_bounds = chop_helices(layers[start:end]) 
        new_start = start+new_bounds[0]
        new_end = start + new_bounds[1] 
        
        tm_helices.append((new_start,new_end))
        #tm_helices.append((start,end))
  else:
    tm_helices = sstruc_bounds
  
  
  
  return tm_helices #, sstruc_bounds



def cylinder_fit_c(n_linepoints, n_unitvector, fragment, num_atoms,path):
  #make the c command - this has the funciton name, then the number of atoms
  # then the linepoints (x y z) and the vector(x y z) from the fit_line and 
  #the helix coordinates (x1 y1 z1 x2 y2 z2 ...) for however many atoms we have 
  
  cylinder = np.zeros([9]) # object for results
  
  command_string = [] #start the list for the command call
  #command_string.append("."+os.sep+"cylinder") #function name
  command_string.append(path+os.sep+"cylinder")
  command_string.append(str(num_atoms)) #number of atoms in fit
  for k in xrange(3): 
    command_string.append(str(n_linepoints[k])) #append the starting point
  for k in xrange(3):
    command_string.append(str(n_unitvector[k])) #append the starting vector
    
  for k in xrange(num_atoms): #for each atom
    command_string.append(str(fragment.get_coords()[k][0])) #append xi
    command_string.append(str(fragment.get_coords()[k][1])) #append yi
    command_string.append(str(fragment.get_coords()[k][2])) #append zi
  
  #run the command string, putting the stdout (which is the linepoints then 
  # the vector, followed by the rmsd, ) into p - a string
  #print command_string
  p = sub.check_output(command_string)
  #p in format %f %f %f %f %f %f %f, so need to parse this
  k=0 # counter for position in p
  
  if p[0] == '*': #an error in the fitting function results in an output of *
    #state where this is, but do not raise warning or exception, as this is one
    # of many possible fits to this helix region
    #print 'fit failed for vector, reverting to least squares fit'
    
    cylinder[3:6] = n_unitvector #return to least squares fit, with rmsds = 1000 and radius = 1000
    cylinder[0:3] = n_linepoints[0]
    cylinder[6] = 1000
    cylinder[7] = 1000
    cylinder[8] = 1000
  
  else: #fit worked
    outputs = p.split(' ') #split the output into individual arguments
    for j in xrange(0,9): #for each of the outputs
      cylinder[j] = float(outputs[j])
    
    #check that vector is going in the same direction as the initial estimate
    if (np.dot(cylinder[3:6], n_unitvector) < 0): # the direction of the fit has been reversed
      cylinder[3:6] = -cylinder[3:6] #switch the direction back to that of the least-squares fit
      
  #take the squared deviations, and return rmsds: 7 is the end result, 8 the 
  # rmsd of the initial guess; 6 is the cylinder radius
  cylinder[7:9] = np.sqrt(cylinder[7:9]/num_atoms)  
  return cylinder 

def start_pymol_script(pdbfile, pdbcode,chain, start, end,  pymol_file_dir = '', 
                       structure_filename = 'none'):
  if structure_filename == 'none':
    print 'error: no filename for structure to be used in pymol script'
    return
  #open a file to write to
  f = open(pymol_file_dir+pdbcode+str(start)+'.py', 'w')
  #instruction to load the pdbfile
  f.write("if '%s' not in cmd.get_names() and '%s' not in cmd.get_names():\n" % (pdbcode, pdbcode[0:4]))
  f.write('\tcmd.load("'+pdbfile+'")\n')
  pdbcode1 = pdbfile.split(os.sep)[-1].split('.')[0]
  #hide what we automatically get
  f.write('\tcmd.hide("everything","'+pdbcode1+'")\n')
  #select the backbond of our helix
  f.write('cmd.select("helix", "resi '+str(start) + ':' + str(end) + ' & name c+ca+n+o & chain ' + chain+'")\n')
  #show the helix as sticks
  f.write('cmd.show("sticks","helix")\n')
  #colours of helix, kink and background
  f.write('cmd.color("green","helix")\n')
  f.write('cmd.bg_color("white")\n')
  f.close()

def write_pymol_kink(n_cylinder, c_cylinder, n_fragment, c_fragment,
                        pdbfile, pdbcode, chain, start, end, kink_pos, kink_ang, pymol_file_dir = ''):
  #open a file to append(!) to
  f = open(pymol_file_dir+pdbcode+str(start)+'.py', 'a')  
  #select the kink residue
  f.write('\n')
  f.write('cmd.select("kink'+str(kink_pos) + '", "resi '+str(kink_pos) + ' & chain ' +chain + '")\n')
  f.write('cmd.color("black","kink'+str(kink_pos)+'")\n')
  #expects cylinder to be a n by 9 array, each row with two points and a radius
  for_pymol_n =  np.zeros([1,7])
  for_pymol_n[0][0:3] = (n_cylinder[0:3] + np.dot(n_cylinder[3:6], \
                                            (n_fragment.get_coords()[0] - \
                                                n_cylinder[0:3])) * n_cylinder[3:6])
  ##print for_pymol_n[0][3:6]
  for_pymol_n[0][3:6] = n_cylinder[0:3] + \
            np.dot(n_cylinder[3:6], \
                  (n_fragment.get_coords()[-1] - \
                  n_cylinder[0:3])) * n_cylinder[3:6]
  for_pymol_n[0][6] = n_cylinder[6]
  
  for_pymol_c =  np.zeros([1,9])
  for_pymol_c[0][0:3] = (c_cylinder[0:3] + np.dot(c_cylinder[3:6], \
                                               (c_fragment.get_coords()[0] - \
                                                c_cylinder[0:3])) * c_cylinder[3:6])
  for_pymol_c[0][3:6] = c_cylinder[0:3] + \
            np.dot(c_cylinder[3:6], \
                  (c_fragment.get_coords()[-1] - \
                  c_cylinder[0:3])) * c_cylinder[3:6]
  for_pymol_c[0][6] = c_cylinder[6]
  
  #make the two cylinders
  f.write('cmd.load_cgo([9.0,')
  for j in xrange(len(for_pymol_n[0])):
    f.write(str(for_pymol_n[0][j])+', ')
  f.write('1,0,0,1,0,0], "n_fit'+str(kink_pos)+'")\n')
  f.write('cmd.set("cgo_transparency", 0.5, "n_fit'+str(kink_pos)+'")\n') 
  f.write('cmd.load_cgo([9.0,')
  for j in xrange(len(for_pymol_c[0])):
    f.write(str(for_pymol_c[0][j])+', ')
  f.write('1,0,0,1,0,0], "c_fit'+str(kink_pos)+'")\n')  
  f.write('cmd.set("cgo_transparency", 0.5, "c_fit'+str(kink_pos)+'")\n') 
  
  #label the kink position
  f.write('''cmd.label("kink'''+str(kink_pos)+''' & name ca", "'Angle = ''' + str(kink_ang)[0:5] + ''' deg '")\n''')
  #and the n cylinder
  f.write('''cmd.label(" i. ''' + str(kink_pos -4) + '''& name c", "'Radius = ''' + str(n_cylinder[6])[0:4] + ''' A '")\n''')
  f.write('''cmd.label(" i. ''' + str(kink_pos -4) + '''& name ca", "'RMSD = ''' + str(n_cylinder[7])[0:4] + ''' A '")\n''')
  f.write('''cmd.label(" i. ''' + str(kink_pos -4) + '''& name n", "' ''' + str(len(n_fragment)) + ''' resi  '")\n''')
  # and the c cylinder
  f.write('''cmd.label(" i. ''' + str(kink_pos +4) + '''& name c", "'Radius = ''' + str(c_cylinder[6])[0:4] + ''' A '")\n''')
  f.write('''cmd.label(" i. ''' + str(kink_pos +4) + '''& name ca", "'RMSD = ''' + str(c_cylinder[7])[0:4] + ''' A '")\n''')
  f.write('''cmd.label(" i. ''' + str(kink_pos +4) + '''& name n", "' ''' + str(len(c_fragment)) + ''' resi '")\n''')  
  f.write('cmd.set("label_position","(3,3,3)")\n')
  
  f.close()     
        
  return

def end_pymol_file(pdbcode, start,pymol_file_dir = ''):
  if os.path.exists(pymol_file_dir+pdbcode+str(start)+'.py'):
    f = open(pymol_file_dir+pdbcode+str(start)+'.py', 'a')  
    #orient on our helix
    
    f.write('cmd.orient("helix")\n')
    f.write('cmd.turn("z",90)\n')
    f.write('cmd.zoom("helix", state=1, buffer =2.0)\n')
    f.write('cmd.ray()\n')
  return

def find_py_index(pdb, pdbpos):
  for i in xrange(len(pdb)):
    if pdb[i].CA.ires == pdbpos:
      break
  return i


def rmsd(crds1, crds2):
  """Returns RMSD between 2 sets of [nx3] np array"""
  assert(crds1.shape[1] == 3)
  assert(crds1.shape == crds2.shape)
  n_vec = np.shape(crds1)[0]
  correlation_matrix = np.dot(np.transpose(crds1), crds2)
  #print correlation_matrix
  v, s, w_tr = np.linalg.svd(correlation_matrix)
  is_reflection = (np.linalg.det(v) * np.linalg.det(w_tr)) < 0.0
  if is_reflection:
    s[-1] = - s[-1]
  E0 = sum(sum(crds1 * crds1)) + \
       sum(sum(crds2 * crds2))
  rmsd_sq = (E0 - 2.0*sum(s)) / float(n_vec)
  rmsd_sq = max([rmsd_sq, 0.0])
  return np.sqrt(rmsd_sq)

def rmsdsum_to_error(rmsdsum):
  if rmsdsum > 0:
    return max(0, np.log(max(rmsdsum - 0.294,0.001))*6.349 +13.15)
  else:
    return 0

####definitions done


def kink_finder(directory, pdb_extension, tem_extension, filename, output_path,
                soluble, display, break_angle, pymol_file_directory, 
                in_out, max_loop_length, user_helices,path):
  
  
  if user_helices != 'none':
    pdb_helices = False
  else:
    pdb_helices = True
  
  find_helices_from_tem = False # this allows the user to input a .tem file 
  #(produced by JOY), and Kink Finder to identify helices from this .tem file
  
  
  #in_out = 'outside' #are we going for the kink to be annotated on the inside or outside?
  in_out = 'inside'
    
  if pymol_file_directory == 'none': #Decide if we are going to write a pymol file
    pymol = False
  else:
    pymol = True
    if not os.path.exists(pymol_file_directory):
      os.makedirs(pymol_file_directory)
  
  num_atoms = 24 #
  helix_vector_length = int(math.ceil(num_atoms / 4))  
  
  #check binary is there
  if not os.path.exists(path+os.sep+'cylinder'):
    print "Cylinder binary not in Kink Finder's directory. See readme.txt for further instructions"
    exit()
  
  #check that the cylinder binary works:
  command_string = [path+os.sep+'cylinder','6',
                     '0','0','0',
                     '0','0','1',
                     '-1','0','0',
                     '0','1','1',
                     '1','0','2',
                     '0','-1','3',
                     '-1','0','4',
                     '0','1','5']
  
  try:
    p = sub.check_output(command_string)
  except:
    print 'Cylinder fitting binary not functioning properly. See readme.txt for further instructions'
    exit()
  
  
  #pdbfiles # we can feed it all of the pdbfiles in the folder
  if filename == 'all':
    pdbfiles = glob("%s%s*%s" % (directory, os.sep,  pdb_extension))
    pdbfiles=sorted(pdbfiles)
  
  else:
    pdbfiles = [filename]
  
  #make a folder for the results
  if not os.path.exists(output_path):
    os.makedirs(output_path)
  
 
  ###############################
  # Open some files to write output to
  ##open a file for writing the all angles to
  angleWriter = csv.writer(open(output_path + 'angles.csv','w'))
  errorWriter = csv.writer(open(output_path + 'errors.csv','w'))
  kinkWriter = csv.writer(open(output_path + 'kinks.csv','w'))
  helixWriter = csv.writer(open(output_path + 'helices.csv','w'))
  
  
  kinkWriter.writerow(["pdb_code","Helix_Start", "Helix_End", "Kink_Position",
          "Kink_Start", "Kink_End", "Kink_Angle", "sequence", "n_radius", 
          "n_rmsd", "c_radius","c_rmsd", 'Estimated_Error'])
  helixWriter.writerow(["pdb_code","Helix_Start", "Helix_End", "Kink_Position",
           "Kink_Angle", "sequence"])
  
  number_formatter = format('{:.3f}')
  display_formatter = format('{:10.3f}')
  res_number_formatter = format('{:6n}')
  
  #for each PDB file
  for pdbfile in pdbfiles:
    #print pdbfile
    pdb_code = pdbfile.split(os.sep)[-1].split('.')[0]
    try:
      assert os.path.exists(pdbfile)
    except:
      print 'File % could not be found' % pdbfile 
      exit() 
    
    #check_for_multiple chains:
    pdb_backbone_atoms = Pdb(pdbfile).get_backbone()
    chains = pdb_backbone_atoms.get_chain_codes()
    
    for chain in chains:
      pdb_code = pdb_code[0:4] + chain
      print 'Analysing chain %s of %s' % (chain, pdb_code)
      pdb = ResidueList(Pdb(pdbfile).get_chain(chain).get_backbone())
      
      
      if find_helices_from_tem:
        temfile = pdbfile[:-len(pdb_extension)] + tem_extension
        try:
          assert os.path.exists(temfile)
        except:
          print 'no',temfile 
          exit()
        tem = Ali(temfile)
        try:
          assert len(pdb) == len(tem[0][0].seq)
        except(AssertionError):
          print 'pdb and tem are different lengths', len(pdb), len(tem[0][0].seq), tem[0][0].seq
          exit()
          
          print len(pdb), tem[0][0].seq
          continue
        
        # Extract secondary structure and membrane layer annotation from TEM file
        sequence = tem[0][0].seq
        sstruc = tem[0]["secondary structure and phi angle"].seq
        #print sstruc
        tm_helices = find_helices(pdb,sstruc,sequence,'', soluble, max_loop_length)
      
      else: #user defined helices
        if user_helices != 'none':
          helix_string = user_helices.split(' ') #parse the string input by users
          helices = []
          for helix in helix_string:
            
            #print helix,
            pdbstart=int(helix.split('-')[0])
            pdbend = int(helix.split('-')[1]) # find where the residue with pbdnumber pdbend is in the list of residues
            #print pdbstart,pdbend, '##',      
            helices.append([pdbstart,pdbend])
          
        else:
          helices = parse_file_for_helices(pdbfile, chain)
          if len(helices) == 0:
            print 'No helix definitions for chain %s in file %s . Either include header in PDB file, or manually specify helix limits' % (chains[0], pdbfile)
            continue
        tm_helices=[]
        
        
        #now covert from the PDB index, to the index in the pdb object
        for helix in helices:
          for j in xrange(len(pdb)):
            if pdb[j].CA.ires == helix[0]:
              start = j
            if pdb[j].CA.ires == helix[1]:
              end = j+1     
          
          tm_helices.append([start,end])
        
        sequence = pdb.get_seq()
        sstruc = 'H' * len(sequence)
      
      
      #============================================================================
      #Loop over helices 
      #============================================================================
      if display:
        print '                   Largest    Largest       '
        print '       First   Last   Kink       Kink     Helix'
        print 'Chain   Resi   Resi   Resi      Angle     Sequence'
      #print tm_helices
      
      for start, end in tm_helices:
        # stop if it is too short - i.e we cannot fit two consecutive cylinders to it
        if end-start < 2*helix_vector_length:
          continue
        #get the coordinates of the helix
        helix=pdb[start:end]
        #check we have all backbone atoms
        if len(helix.get_coords())!= 4*(end-start):
          print "Only %i of %i backbone atoms in coordinate file for helix starting at residue %i of %s. Kink Finder needs all backbone atoms (CA, C, O, N)" % ( len(helix.get_coords()),4*(end-start),start, pdb_code)
          continue 
        
        
        #lets initialise some things
        maxangle=0.0 #the biggest angle
        maxpos=-1 #position of the kinks
        angles = [] #angles of the helix
        angles2 = [] #wobble angles
        #initialize some arrays for unitvectors and cylinder results
        unitvector= np.zeros([end-start-helix_vector_length+1,3]) # an array of vectors
        cylinder = np.zeros([end-start-helix_vector_length+1,9]) # an array of fitted cylinders
        
        #==========================================================================
        # loop over the helix, calculating initial vectors
        #==========================================================================
        for i in xrange(end-start-helix_vector_length+1):
          #choose our fragment to fit
          fragment = helix[i:i+helix_vector_length]
          #simple least squares to get a sensible starting point       
          n_unitvector, n_linepoints = fit_line(fragment.get_coords()[0:21])
          #now cylinder fit, using the least squares 
          cylinder[i] = cylinder_fit_c(n_linepoints[0], n_unitvector,  fragment, num_atoms,path)
          
          #========================================================================
          # Deal with situations where the initial cylinder fit is poor. This 
          # occurs where the helix is distorted, and so a least squares fit to 21 
          # atoms does not give a good approximation, so the cylinder fitting routine
          # gets stuck in a local minimum. Pi and tight turns are better approximated 
          # by 16 (tight), 26(pi) and 28 (tight) atom fits. Using a longer region to
          # fit to also gives a better approximation of the helix axis, so a 36 atom
          # least-squares fit is also included gives
          # a 9 residue section, whic 
          #========================================================================
          
          #if the fit is not great, try a series of other starting points - this uses
          #26, 16, 28 and 36 atoms to fit a least squares line to. Then uses that for 
          #a starting point for the 24 atom cylinder fit
          if cylinder[i][7] > 0.31:
            #make a new object, to include the cylinder fits we do
            cylinder2 = np.ones([6,9])*10
            #try a 26 atom starting fit, in case of pi helix
            f_start = max(0,(i-1))
            f_end = f_start + 7
            fragment2 = helix[f_start:f_end]
            n_unitvector, n_linepoints = fit_line(fragment2.get_coords()[0:26])
            cylinder2[0] = cylinder_fit_c(n_linepoints[0], n_unitvector,  fragment, 24,path)
            
            # try a 16 atom starting fit, in case of 3_10 
            f_start = i+1
            f_end = i+helix_vector_length-1
            fragment2 = helix[f_start:f_end]
            n_unitvector, n_linepoints = fit_line(fragment2.get_coords()[0:16])
            cylinder2[1] = cylinder_fit_c(n_linepoints[0], n_unitvector,  fragment, 24,path)
            
            #try a 30 atom fit
            if i == 0:
              f_start = 0
            else:
              f_start = min((i-1),(end-start-helix_vector_length-2))
            f_end = f_start+8
            fragment2 = helix[f_start:f_end]
            n_unitvector, n_linepoints = fit_line(fragment2.get_coords()[0:28])
            cylinder2[2] = cylinder_fit_c(n_linepoints[0], n_unitvector,  fragment, 24,path)
            
            #try a 36 atom fit
            if i == 0:
              f_start = 0
            else:
              f_start = min((i-1),(end-start-helix_vector_length-3))
            f_end = f_start+9
            fragment2 = helix[f_start:f_end]
            n_unitvector, n_linepoints = fit_line(fragment2.get_coords()[0:36])
            cylinder2[3] = cylinder_fit_c(n_linepoints[0], n_unitvector,  fragment, num_atoms,path)
            
            #try the fitted axis from the previous section of the helix 
            #(if this is not the first section)
            if i !=0:
              cylinder2[4] = cylinder_fit_c(cylinder[i-1][0:3],cylinder[i-1][3:6],fragment,num_atoms,path)
            
            # and lets include the initial fit in this, in case it is the best
            cylinder2[5] = cylinder[i]
            
            #======================================================================
            # now want the fit with the best rmsd, given that it has a sensible diameter
            #======================================================================
            okay_radius = (abs(cylinder2[:,6]-2.0)<0.3) # we are within the range of good radii 
            okay_rmsd = cylinder2[:,7]<0.38 # we have a low rmsd
            okay = ((okay_radius + okay_rmsd) == 2)
            cylinder_rms_okay = cylinder2[okay]
            
            if len(cylinder_rms_okay) != 0:
              #if we have something(s) that are within both of these constraints, 
              #pick the one with the best (lowest) rmsd
              best_cylinder2 = cylinder_rms_okay[cylinder_rms_okay[:,7].argmin()]
              #now replace the original cylinder fit with this one
              cylinder[i]=best_cylinder2 
              #print best_cylinder2
            else:
              #now we have no fit with a sensible radius, so pick the one with the best score
              #the score is |radius -2 | + 5*(rmsd - 0.25)
              scores = abs(cylinder2[:,6]-2.0) + 5*(cylinder2[:,7]-0.25) 
              cylinder[i] = cylinder2[scores.argmin()]
          
          #========================================================================
          # We now have a cylinder axis for each sliding window of 6 residues
          #========================================================================
          
          
          unitvector[i] = cylinder[i][3:6] #transfer to unitvector
        
        #==========================================================================
        # now run back and forth down the helix, trying the next fit as a starting point
        # do this until it has been done 10 times, or there has been no change in the fits
        #==========================================================================
        return_changes =1
        iterations = 0
        while return_changes > 0 and iterations < 10:
          iterations +=1
          return_changes = 0
          #print return_changes, 'hi'
          #now go back and try the next fit for each
          for i in xrange((end-start-helix_vector_length+1-2),-1,-1):
            #choose our fragment to fit
            fragment = helix[i:i+helix_vector_length]
            #fit using the fit from the previous chunk
            cylinder3 = cylinder_fit_c(cylinder[i+1][0:3],cylinder[i+1][3:6],fragment,num_atoms,path)
            
            #if the r is sensible, and the rmsd is better by a significant margin,
            if (abs(cylinder3[6]-2.0) < 0.3) and cylinder[i][7] - cylinder3[7] > 0.01:
              #print return_changes, cylinder3[6:8], cylinder[i][6:8]
              cylinder[i] = cylinder3
              return_changes +=1
              
            elif (abs(cylinder[i][6]-2.0) > 0.3) or cylinder[i][7] > 0.38: 
              #else if we have a really quite bad previous fit
              score_old = abs(cylinder[i][6]-2.0) + 5*(cylinder[i][7]-0.25)
              score_new = scores = abs(cylinder3[6]-2.0) + 5*(cylinder3[7]-0.25)
              if score_new < score_old:
                cylinder[i] = cylinder3
                return_changes +=1
          
          #and try the following fit:
          for i in xrange(1,end-start-helix_vector_length+1):
            #choose our fragment to fit
            fragment = helix[i:i+helix_vector_length]
            #fit using the fit from the previous chunk
            cylinder3 = cylinder_fit_c(cylinder[i-1][0:3],cylinder[i-1][3:6],fragment,num_atoms,path)
            
            #if the r is sensible, and the rmsd is better by a significant margin,
            if (abs(cylinder3[6]-2.0) < 0.3) and cylinder[i][7] - cylinder3[7] > 0.01:
              #print return_changes, cylinder3[6:8], cylinder[i][6:8]
              cylinder[i] = cylinder3
              return_changes +=1
              
            elif (abs(cylinder[i][6]-2.0) > 0.3) or cylinder[i][7] > 0.38: 
              #else if we have a really quite bad previous fit
              score_old = abs(cylinder[i][6]-2.0) + 5*(cylinder[i][7]-0.25)
              score_new = scores = abs(cylinder3[6]-2.0) + 5*(cylinder3[7]-0.25)
              if score_new < score_old:
                cylinder[i] = cylinder3
                return_changes +=1
                
        
        #==========================================================================
        # #Now, have best 6 residue fits
        # #Want to see if longer fits may be better
        #==========================================================================
        number_of_possible_angles = len(helix)-2*helix_vector_length+1
        new_cylinder_n = np.zeros([number_of_possible_angles,9]) #for fits on the n side of the kink
        new_cylinder_c = np.zeros([number_of_possible_angles,9]) #for fits on the c side of the kink
        ############### try longer helices everywhere: 3/8/12
        
        #print start, end, len(helix), helix_vector_length, number_of_possible_angles
        for j in xrange(number_of_possible_angles): #for each possible kink point
          #initialise an array
          #want to go from 5 to 10, but only if that is possible - so take min of 6
          # and j, which is 0 if at the first point in the helix
          n_cylinders = np.zeros([min(6,j+1),9]) #so for up to 6 different fits
          #n_cylinders[1]
          
          for i in xrange(min(6,j+1)): # +1 due to python counting (if 
            # it was 0, then there would be no calculation done.
            #work out some longer vectors - 
            #take the fragment
            n_fragment = helix[(j-i):(j+helix_vector_length)]  
            
            #fit, using the original vector as a start
            #n_linepoints, n_unitvector,  fragment, num_atoms
            number_of_fragment_atoms = (i+6)*4
            n_cylinders[i]=cylinder_fit_c(cylinder[j][0:3], 
                                          cylinder[j][3:6],
                                          n_fragment, number_of_fragment_atoms,path)
            #print n_cylinders[i]
          #print n_cylinders
          new_cylinder_n[j] = n_cylinders[n_cylinders[:,7].argmin()] #best is one with smallest rmsd
          
          c_cylinders = np.zeros([min(6,number_of_possible_angles-j),9]) 
          #print i, 'c_cylinders', (len(helix)-j-1)
          for i in xrange(min(6,(number_of_possible_angles-j))):
            #work out some longer vectors
            #take the fragment
            #print i, 'c_cylinders', (len(helix)-j)
            c_fragment = helix[(j+helix_vector_length):
                              (j+2*helix_vector_length+i)]
                
            #fit, using the original vector as a start
            #c_linepoints, n_unitvector,  fragment, num_atoms
            number_of_fragment_atoms = ((i+6)*4)
            
            
            c_cylinders[i]=cylinder_fit_c(cylinder[j+helix_vector_length][0:3], 
                                          cylinder[j+helix_vector_length][3:6],
                                          c_fragment, number_of_fragment_atoms,path)
          #print c_cylinders
          #print n_cylinders
          try:
            new_cylinder_c[j] = c_cylinders[c_cylinders[:,7].argmin()]
          except(ValueError):
            print len(helix), j, number_of_possible_angles,c_cylinders,min(6,(number_of_possible_angles-j))
            print 'error! Should exit'
            #best is one with smallest rmsd
        
        
        
        #==========================================================================
        # We now have the cylinder fits
        # Now, Calculate a kink angle for each residue
        #==========================================================================
        # initialise:
        wobble_angles = []
        outside_helix =[] 
        angles = []
        angle_error_list = []
        
        for i in xrange(len(new_cylinder_c)):
          angles.append(abs(math.degrees(angle(new_cylinder_n[i][3:6], 
                                     new_cylinder_c[i][3:6]))))
          #calculate error estimate from sum of rmsds
          angle_error_list.append(rmsdsum_to_error(new_cylinder_c[i][7] + new_cylinder_n[i][7]))
          #and the 2nd (wobble) angle is given by: 
          wobble_angles.append(wobble2(helix[i+helix_vector_length-1].CA.xyz,
                                           new_cylinder_n[i][0:3],
                                           new_cylinder_n[i][3:6],
                                           new_cylinder_c[i][3:6]))
          angles2.append(wobble_angles[i])
          #near to 0 for outside the kink
          outside_helix.append(abs(wobble_angles[i]-180))
          if maxangle < angles[i]:
            maxangle = angles[i]
            max_position_helix = i
        maxpos = max_position_helix+helix_vector_length-1+start
        ## want to cut the helix if it contains a kink angle > 100
        if maxangle > break_angle:
          # make 2 new helices - start: max pos, maxpos:end
          tm_helices.append([start, maxpos-1])
          tm_helices.append([maxpos+1, end])
          #continue to next helix, without including this in the written documents
          continue 
        #===========================================================================
        # this is where we identify the kink point(s) 
        #===========================================================================
        ##list the angles, and sort them    
        sorted_angles = list(angles) # make a copy of the list
        sorted_angles.sort() # sort it
        sorted_angles.reverse() # biggest first
        
        helix_kink_angle = sorted_angles[0]
        #work out the position of the angles in decreasing order
        kink_indices = []
        for i in xrange(len(sorted_angles)):
          kink_indices.append(angles.index(sorted_angles[i]))
          #this is a list of the index of the largest to smallest kink angle
        
        #now go down this index list picking out kinks
        kink_pos = list([kink_indices[0]]) #the first entry is obviously a kink
        k=1 #counter
        #step down the index list
        for i in xrange(len(kink_indices)-1):
          #stop if the angle is less than 10 # perhaps these should be a while loop,
          #but it didnt seem to work
          if angles[kink_indices[k]] > 10.0 :
            #have to be at least n away from the last kink
            #calculate how close this residue is to an exsisting kink
            min_distance = 1000
            for j in xrange(len(kink_pos)):
              #calculate distance from each other kink
              distance = abs(kink_indices[k] - kink_pos[j])
              #take the minimum distance
              min_distance = min(min_distance,distance)
            #providing we are not within 6 of an existing kink     
            if min_distance > 6:
              # and there is an angle under 10 degrees in between
              straight_in_between = 1
              for j in xrange(len(kink_pos)):
                if min(angles[min(kink_indices[j],
                                  kink_indices[k]):
                              max(kink_indices[j],
                                  kink_indices[k])]) > 10.0:
                  straight_in_between = 0 # there is no angle between this proposed kink and any of the others which is under 10 degrees
                #then we have another kink!
              if straight_in_between:
                kink_pos.append(kink_indices[k])
            #add one to our counter
          k=k+1
        
        
        ##start writing a pymol file 
        if pymol:
          start_pymol_script(pdbfile, pdb_code,chain, pdb[start].CA.ires, 
                       pdb[end-1].CA.ires, pymol_file_dir=pymol_file_directory,
                       structure_filename = pdbfile)
        
        # for each kink
        #===========================================================================
        #  work out some longer vectors, and compare the rmsds - picking the fit with
        # the smallest rmsd
        # the indices for the following things are as follows:
        #
        #  RRRRRRRRRRRRRRRRRRR
        #  0123456789...      Residues, annotatation
        #  000000         }
        #   111111        }
        #    222222       }
        #     333333      }cylinders
        #      ....       }
        #        666666   }
        #           ..... }
        #       0123456789...  Angles
        #  ..0123456789....  C-alpha
        #===========================================================================    
        #print start, end, kink_pos
        a = 0
        for maxang_helix in kink_pos:  #maxang_helix indices are offset from the helix indices by 6 
          #initialise an array
          best_n_cylinder = new_cylinder_n[maxang_helix]
          best_c_cylinder = new_cylinder_c[maxang_helix]
          
          n_fragment = helix[maxang_helix:(maxang_helix+helix_vector_length)]
          c_fragment = helix[(maxang_helix+helix_vector_length):(maxang_helix+2*helix_vector_length)]
          
          #print kink_indices
          #print len(helix), (maxang_helix-helix_vector_length),maxang_helix,(maxang_helix+helix_vector_length)
          #print maxang_helix, maxang_helix+2*helix_vector_length, len(helix)
          #print helix[0]
          #print n_fragment.get_coords()
          #print c_fragment.get_coords()
          
          ################
          # Work out the kink_angle
          kink_angle = math.degrees(angle
                                    (best_n_cylinder[3:6], best_c_cylinder[3:6]))
          #work out some wobble angles (+- six residues)
          wobble_angles_kink = np.zeros([12])
          for i in xrange(helix_vector_length):
            if (i+maxang_helix)<0 :
              wobble_angles_kink[i] = 400
            else:
              wobble_angles_kink[i] = (wobble2(helix[i+maxang_helix].CA.xyz,
                                        best_n_cylinder[0:3],
                                        best_n_cylinder[3:6],
                                        best_c_cylinder[3:6]))
          for i in xrange(helix_vector_length,2*helix_vector_length):
            #these must have the c and v  vectors reversed. As now effectively looking 
            #from the other end of the helix, these are 360-angle
            
            if len(helix)<= (i+maxang_helix):
              wobble_angles_kink[i] = np.nan
            else:
              wobble_angles_kink[i] = (360-wobble2(helix[i+maxang_helix].CA.xyz,
                                        best_c_cylinder[0:3],
                                        -1*best_c_cylinder[3:6],
                                        -1*best_n_cylinder[3:6]))
          
          
          #position_of_kink_in_section = min(6,maxang_helix)
          
          ## outside is wobble = 180 
          #look at the wobble angles around the kink, and pick the one nearest 180
    #      position_correction = ((abs(180 - wobble_angles_kink
    #                                           [(helix_vector_length - 2):
    #                                            (helix_vector_length + 2)])).argmin())
          #=========================================================================
          # This is repositioning the kink
          #=========================================================================
          if in_out == 'outside':
            position_correction = np.nanargmin(abs(180 - wobble_angles_kink
                                               [(helix_vector_length - 2):
                                                (helix_vector_length + 2)]))
          elif in_out == 'inside':
            position_correction = np.nanargmax(abs(180 - wobble_angles_kink
                                               [(helix_vector_length - 2):
                                                (helix_vector_length + 2)]))
          
          corrected_kink_position_protein = maxang_helix+helix_vector_length-1+start-1+position_correction
          uncorrected_kink_position_protein = maxang_helix+helix_vector_length-1+start
          
          
          if a ==0:
            helix_corrected_biggest_kink_position_protein = corrected_kink_position_protein
            a+=1
          
          #========================================================================
          # Write the information to the Kink file
          #========================================================================
          
          kink_start = max(corrected_kink_position_protein - helix_vector_length,0)
          kink_end = min(kink_start+2*helix_vector_length+1,len(pdb))
          #write kink, starting with general info about the kink, then, specifically about the kink
          
          
          kinkWriter.writerow([pdb_code, pdb[start].CA.ires, pdb[end-1].CA.ires, 
                pdb[corrected_kink_position_protein].CA.ires, pdb[kink_start].CA.ires, 
                pdb[kink_end-1].CA.ires,number_formatter.format(kink_angle),
                sequence[kink_start:kink_end],
                number_formatter.format(best_n_cylinder[6]),
                number_formatter.format(best_n_cylinder[7]), 
                number_formatter.format(best_c_cylinder[6]), 
                number_formatter.format(best_c_cylinder[7]),
                number_formatter.format(rmsdsum_to_error(best_n_cylinder[7] + best_c_cylinder[7]))])
          
          #========================================================================
          # Write the kink information to the pymol file
          #========================================================================      
          
          if pymol:
            write_pymol_kink(best_n_cylinder, best_c_cylinder, n_fragment, 
                           c_fragment, pdbfile, pdb_code, chain, pdb[start].CA.ires, 
                           pdb[end-1].CA.ires, 
                           pdb[corrected_kink_position_protein].CA.ires, kink_angle,
                           pymol_file_dir=pymol_file_directory)
          
        
        #==========================================================================
        # Now, print the helix info, displaying it if desired
        #==========================================================================
        
        #print the info about the helix
        if display == True:
          
          #print pdb_code, maxangle, pdb[maxpos].CA.ires, angles
          print pdb_code, res_number_formatter.format(pdb[start].CA.ires), \
              res_number_formatter.format(pdb[end-1].CA.ires), \
              res_number_formatter.format(pdb[helix_corrected_biggest_kink_position_protein].CA.ires), \
              display_formatter.format(helix_kink_angle), '   ', \
              sequence[start:end]
        
        
        helixWriter.writerow([pdb_code, pdb[start].CA.ires, pdb[end-1].CA.ires, 
              pdb[helix_corrected_biggest_kink_position_protein].CA.ires, 
              number_formatter.format(helix_kink_angle),
              sequence[start:end]])
        
        angle_strings = [pdb_code+str(pdb[start].CA.ires),'0','0','0','0','0'] #the code of the helix + first residue number, then 0 for the first 5 residues where no angle has been calculated.
        error_strings = [pdb_code+str(pdb[start].CA.ires),'0','0','0','0','0']
        for angle1 in angles:
          angle_strings.append(number_formatter.format(angle1))
        for error in angle_error_list:
          error_strings.append(number_formatter.format(error))
        angle_strings.extend(['0','0','0','0','0','0']) #no angle for last 6 residues
        error_strings.extend(['0','0','0','0','0','0']) #no angle for last 6 residues
        
        angleWriter.writerow(angle_strings)
        errorWriter.writerow(error_strings)
        #==========================================================================
        # Close the pymol file for this helix
        #==========================================================================
        end_pymol_file(pdb_code, pdb[start].CA.ires,pymol_file_dir=pymol_file_directory)
        
        
  return

def parse_file_for_helices(filename, chain_code):
  file_pdb = open(filename,'r')
  helices = []
  for line in file_pdb:
    if (line[0:5] == 'HELIX') and line[39] == '1' and line[19] == chain_code: #an alpha helix, in the correct chain
        helices.append([int(line[21:25]),int(line[33:37])])
  
  return helices
