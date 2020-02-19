#!/usr/bin/python

import sys
import os
import math

import numpy

from prosci.common import join, average
from prosci.util.ali import Ali
from prosci.util.pdb import Pdb
from prosci.util.pdb3d import fit_line, angle
from prosci.util.gaps import isGap, count_gaps, deGappify


# Annotation value that signifies "in the middle layer of the membrane"
def isTransMembrane(char):
  return char == "T"

def findTMregion(seq, minsize):
  slices = []
  
  # Search for TM regions longer than minsize
  start=-1
  for i,v in enumerate(seq):
    if isGap(v):
      continue
    
    if isTransMembrane(v):
      if start < 0:
        start = i
    else:
      if start >= 0:
        if i - start >= minsize:
          slices.append([start, i])
        start = -1
  
  # Special case where TM region is at end of sequence
  if start >= 0:
    i = len(seq)
    while isGap(seq[i-1]) and i>start:
      i -= 1
    
    assert isTransMembrane(seq[i-1])
    if i - start >= minsize:
      slices.append([start, i])
  
  return slices


def extendTMregion(slices, seq):
  nativeslices = [x[:] for x in slices]
  nativeslices.insert(0, [None, 0])
  nativeslices.append([len(seq), None])
  
  for sl_num, sl in enumerate(slices):
    start, end = sl
    left_limit = nativeslices[sl_num][1]
    right_limit = nativeslices[sl_num+2][0]
    
    # extend to the left
    i = start-1
    while i >= left_limit:
      if isTransMembrane(seq[i]):
        start = i
      elif not isGap(seq[i]):
        break
      i-=1
    
    # extend to the right
    i = end
    while i < right_limit:
      if isTransMembrane(seq[i]):
        end = i+1
      elif not isGap(seq[i]):
        break
      i+=1
    
    sl[:] = start, end


def reduceTMregion(slices, seq, minsize):

  sl_num=0
  while sl_num < len(slices):
    start, end = slices[sl_num]
    
    # reduce left end
    i = start
    while i < end:
      if not isTransMembrane(seq[i]):
        start = i+1
      else:
        break
      i+=1
    
    # reduce right end
    i = end-1
    while i >= start:
      if not isTransMembrane(seq[i]):
        end = i
      else:
        break
      i-=1
    
    if end - start < minsize:
      del slices[sl_num]
    else:
      slices[sl_num][:] = start, end
      sl_num+=1


def syncTMregion(bigslices, smallslices):
  def includes(a,b):
    return a[0] <= b[0] and a[1] >= b[1]
  
  deleted=0
  
  i=0
  while i < len(bigslices):
    bs = bigslices[i]
    isok = False
    for ss in smallslices:
      if includes(bs, ss):
        isok = True
        break
    if not isok:
      del bigslices[i]
      deleted += 1
    else:
      i += 1
  
  i=0
  while i < len(smallslices):
    ss = smallslices[i]
    isok = False
    for bs in bigslices:
      if includes(bs, ss):
        isok = True
        break
    if not isok:
      del smallslices[i]
      deleted += 1
    else:
      i += 1
  
  assert len(bigslices) == len(smallslices)
  
  if deleted:
    sys.stderr.write("Deleted %d tm regions.\n"%(deleted))
  

def diffTM(slices, seq1, seq2):
  scorelist=[]
  for start, end in slices:
    diff=0
    for i in xrange(start, end):
      if seq1[i] != seq2[i]:
        diff += 1
    scorelist.append(diff)
  
  if scorelist:
    avgscore = float(sum(scorelist)) / len(scorelist)
  else:
    avgscore = 0.0
  
  return avgscore, scorelist


def columnindex2residueindex(slices, seq):
  resslices = []
  for start, end in slices:
    start -= count_gaps(seq, 0, start)
    end   -= count_gaps(seq, 0, end)
    resslices.append([start, end])
  return resslices


#~ def slice2range(slices):
  #~ indeces=[]
  #~ for start,end in slice:
    #~ indeces.extend(xrange(start, end))
  #~ return indeces




def project_point_onto_vector(vector_anchor, unitvector, point):
  hypothenuse = point - vector_anchor
  Lh = numpy.linalg.norm(hypothenuse)
  
  theta = angle(unitvector, hypothenuse)
  L2 = Lh*math.cos(theta)
  
  projected_point = unitvector * L2 + vector_anchor
  
  #~ print "len(hypothenuse)", Lh
  #~ print "len(cathete)", L2
  #~ print "theta", rad2deg(theta)
  #~ print "angle2", rad2deg(angle(point - projected_point, hypothenuse))
  #~ print "angle3", rad2deg(angle(point - projected_point, unitvector))
  
  return projected_point


#~ def calculate_absolute_angles(model):
  #~ from prosci.util.pdb3d import fit_line, angle
  #~ import math
  #~ import numpy
  #~ 
  #~ zaxis = numpy.array([0.0, 0.0, 1.0])
  #~ model_axis_unitvector,  model_axis_points  = fit_line(model)
  #~ 
  #~ tilt_angle = angle(zaxis, model_axis_unitvector)
  #~ 
  #~ # get coordinates of middle CA atom
  #~ #
  #~ model_middle_atom_coords  = get_middle_CA_coords(model)
  #~ 
  #~ # build a triangle with the fragment axis and the middle CA atom
  #~ #
  #~ model_axis_point  = project_point_onto_vector(model_axis_points[0],  model_axis_unitvector,  model_middle_atom_coords)
  #~ 
  #~ rS = model_middle_atom_coords - model_axis_point
  #~ rSxA = numpy.cross(rS/numpy.linalg.norm(rS), model_axis_unitvector)
  #~ 
  #~ # now project xaxis onto plane formed by rS and rSxA ...


class TMAngleError(ValueError):
  pass


def rad2deg(angle):
  return angle*180.0/math.pi


class TMFragment:
  zaxis = numpy.array([0.0, 0.0, 1.0])
  
  def __init__(self, atoms, update=True):
    if isinstance(atoms, TMFragment):
      self.atoms = atoms.atoms.deep_copy()
    else:
      self.atoms = atoms.deep_copy()
      self.update()
    if update:
      self.update()
  
  
  def update(self):
    self.axis_unitvector, self.axis_points = fit_line(self.atoms)
    
    # get tilt angles
    #
    self.tilt_angle = angle(TMFragment.zaxis, self.axis_unitvector)
    
    # get coordinates of middle CA atom
    #
    self.middle_atom_coords = self.get_middle_CA_coords()
  
    # build a triangle with the fragment axis and the middle CA atom
    #
    self.middle_axis_point = project_point_onto_vector(self.axis_points[0], self.axis_unitvector, self.middle_atom_coords)
    
    # get vector used for calculating the rotation angle
    self.rS = self.middle_atom_coords - self.middle_axis_point
    self.rS_points = numpy.array([ self.middle_axis_point, self.middle_atom_coords])
    
    #~ rotvec = - numpy.array([ self.axis_unitvector[0], self.axis_unitvector[1] ])
    #~ #rotvec /= numpy.linalg.norm(rotvec)
    #~ rotvec2 = numpy.array([ self.rS[0], self.rS[1] ])
    #~ #rotvec2 /= numpy.linalg.norm(rotvec2)
    #~ self.rotation_angle = angle(rotvec2, rotvec)
    
    #rSxA = numpy.cross(rS/numpy.linalg.norm(rS), model_axis_unitvector)
  
  
  def isCorrectOrientation(self):
    temp = abs(self.tilt_angle) % math.pi
    if temp > math.pi/2.0:
      temp = math.pi - temp
    if temp > math.pi/4.0:
      return False
    return True
  
  def get_middle_CA_coords(self):
    # get coordinates of middle CA atom
    #
    CAs = self.atoms.get_CA()
    middle_atom = CAs[len(CAs)/2]
    return numpy.array([middle_atom.x, middle_atom.y, middle_atom.z])
  
  def transform(self, rotation_matrix, rotation_centre, move_to=None):
    from prosci.util.pdb3d import get_coords, set_coords, rotate, translate
    
    if move_to is None:
      move_to = rotation_centre
    
    c = get_coords(self.atoms)
    translate(c, -rotation_centre)
    rotate(c, rotation_matrix)
    translate(c, move_to)
    
    newFrag = TMFragment(self)
    set_coords(newFrag.atoms, c)
    newFrag.update()
    return newFrag
  
  def overlay_onto(self, other):
    from prosci.util.pdb3d import vectors2rotation_matrix
    rotmat = vectors2rotation_matrix(self.axis_unitvector, other.axis_unitvector)
    return self.transform(rotmat, self.middle_axis_point, other.middle_axis_point)

  def get_relative_rotation(self, other):
    from prosci.util.pdb3d import vectors2rotation_matrix
    rotmat = vectors2rotation_matrix(other.axis_unitvector, self.axis_unitvector)
    new_rS = numpy.dot(other.rS, rotmat)
    #print other.rS, new_rS
    return angle(self.rS, new_rS)

  def get_relative_tilt(self, other):
    return self.tilt_angle - other.tilt_angle



def plot_fragments(native_frag, model_frag):
    import matplotlib.pyplot as plt
    import mpl_toolkits.mplot3d as m3d
    from prosci.util.pdb3d import get_coords
    
    native_coords = get_coords(native_frag.atoms.get_CA())
    model_coords = get_coords(model_frag.atoms.get_CA())
    
    #native_hypothenuse = numpy.array([ native_frag.axis_points[0], native_frag.middle_atom_coords ])
    #model_hypothenuse = numpy.array([ model_frag.axis_points[0], model_frag.middle_atom_coords ])
    
    ax = m3d.Axes3D(plt.figure())
    ax.plot3D(*native_coords.T, color="r")
    ax.plot3D(*native_frag.axis_points.T, color="r", linewidth=5)
    ax.plot3D(*native_frag.rS_points.T, color="r", linewidth=3)
    #ax.plot3D(*native_hypothenuse.T, color="r", linewidth=3)
    ax.plot3D(*model_coords.T, color="b")
    ax.plot3D(*model_frag.axis_points.T, color="b", linewidth=5)
    ax.plot3D(*model_frag.rS_points.T, color="b", linewidth=3)
    #ax.plot3D(*model_hypothenuse.T, color="b", linewidth=3)
    ax.set_aspect('equal', 'datalim') # no effect??
    plt.show()



class ModelComparer:
    def __init__(self, alignment, seqalign=True):
      if isinstance(alignment, Ali):
        self.alignment = alignment
      else:
        self.alignment = Ali(alignment)
      
      if seqalign:
        self.alignment.align(degappify=True) # we want alignment by residue equivalence, not by structural similarity
      
      self.structures = {}
    
    
    def load_structure(self, code, fname):
      self.structures[code] = Pdb(code, file(fname))
    
    def load_structures(self):
      for eg in self.alignment:
        fname = self.alignment[nativeid].getMasterEntry().getStructureFilename()
        if fname:
          code = eg.getCode()
          self.structures[code] = Pdb(code, file(fname))
      
    
    #~ def score_tm_shift(self, nativeid, modelid, minlength=7):
      #~ native = self.alignment[nativeid]["membrane layer"]
      #~ model = self.alignment[modelid]["membrane layer"]
      #~ 
      #~ assert len(native.seq) == len(model.seq)
      #~ 
      #~ tm_region = findTMregion(native.seq, minlength)
      #~ extendTMregion(tm_region, model.seq)
      #~ 
      #~ return tm_region, diffTM(tm_region, native.seq, model.seq)
    
    
    def score_tm_segments(self, nativeid, modelid, minlength=7, debug=False):
      native = self.alignment[nativeid]["membrane layer"]
      model = self.alignment[modelid]["membrane layer"]
      
      assert len(native.seq) == len(model.seq)
      
      if not self.structures:
        self.load_structures()
      
      tm_region  = findTMregion(native.seq, minlength)
      tm_region2 = [ x[:] for x in tm_region ]
      reduceTMregion(tm_region,  model.seq, minlength)
      extendTMregion(tm_region2, model.seq)
      syncTMregion(tm_region2, tm_region)
      
      tm_region_native = columnindex2residueindex(tm_region, native.seq)
      tm_region_model  = columnindex2residueindex(tm_region, model.seq)
      
      assert len(tm_region_native) == len(tm_region_model)
      
      tm_fragments_native = [ self.structures[nativeid].get_residue_slice(start, end-start) for start,end in tm_region_native ]
      tm_fragments_model  = [  self.structures[modelid].get_residue_slice(start, end-start) for start,end in tm_region_model  ]
      
      assert len(tm_fragments_native) == len(tm_fragments_model)
      assert len(tm_fragments_native) == len(tm_region)
      
      angles = []
      i=0
      while i < len(tm_fragments_native):
        frag_nat = tm_fragments_native[i]
        frag_mod = tm_fragments_model[i]
        
        alistart, aliend = tm_region[i]
        assert frag_nat.get_seq() == deGappify(self.alignment[nativeid].getMasterEntry().seq[alistart:aliend])
        assert frag_mod.get_seq() == deGappify(self.alignment[modelid].getMasterEntry().seq[alistart:aliend])
        assert frag_nat.rescount() == frag_mod.rescount()
        
        frag_nat = TMFragment(frag_nat)
        if not frag_nat.isCorrectOrientation():
          del tm_region[i]
          del tm_region2[i]
          del tm_region_native[i]
          del tm_region_model[i]
          del tm_fragments_native[i]
          del tm_fragments_model[i]
          continue
        
        frag_mod  = TMFragment(frag_mod)
        tilt_angle = frag_nat.get_relative_tilt(frag_mod)
        rotation_angle = frag_nat.get_relative_rotation(frag_mod)
        angles.append( (tilt_angle, rotation_angle) )
        i += 1

        #frag_mod2  = frag_mod.overlay_onto(frag_nat)
        #print rad2deg(angle(frag_nat.rS, frag_mod2.rS))
        #print rad2deg(angle(frag_nat.rS, frag_mod.rS))
        if debug:
          plot_fragments(frag_nat, frag_mod)

      
      shifts = diffTM(tm_region2, native.seq, model.seq)
      
      return tm_region, angles, tm_region2, shifts



if __name__ == "__main__":

    # filename refers to a .tem file.
    if len(sys.argv) != 6:
      print "usage: python tm.py <tem_file> <native_id> <model_id> <native_structure> <model_structure>"
      sys.exit(1)
    filename, nativeid, modelid, nativestruc, modelstruc = sys.argv[1:]
    comparer = ModelComparer(file(filename))
    
    #~ tm_region_shift, (savg, slist) = comparer.score_tm_shift(nativeid, modelid)
    #~ print "%s %2d %5.2f %s"%(filename, len(slist), savg, join(",", slist))
    
    comparer.load_structure(nativeid, nativestruc)
    comparer.load_structure(modelid, modelstruc)

    # Calculate angles and shifts
    tm_region_angles, angles, tm_region_shifts, (savg, slist) = comparer.score_tm_segments(nativeid, modelid)
    
    # Covert radians to degrees
    for i,a in enumerate(angles):
      angles[i] = [ rad2deg(x) for x in a ]
    
    # Get average tilt and rotation angles
    average_tilt = average([x[0] for x in angles])
    average_rot  = average([x[1] for x in angles])
    
    
    print "%s %d %6.2f %6.2f"%(filename, len(angles), average_tilt, average_rot)
    for i,a in enumerate(angles):
      alistart, aliend = tm_region_angles[i]
      print "\tTM %2d"%(i+1)
      print "\tnative sequence:   %s"%(comparer.alignment[nativeid].getMasterEntry().seq[alistart:aliend])
      print "\tmodel  sequence:   %s"%(comparer.alignment[modelid].getMasterEntry().seq[alistart:aliend])
      print "\tnative layer:      %s"%(comparer.alignment[nativeid]["membrane layer"].seq[alistart:aliend])
      print "\tmodel  layer:      %s"%(comparer.alignment[modelid]["membrane layer"].seq[alistart:aliend])
      print "\trelative tilt:     %6.2f degrees"%(a[0])
      print "\trelative rotation: %6.2f degrees"%(a[1])
      print
    
