#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Residue class represents a single amino acid in a PDB file.
#
# Author: Sebastian Kelm
# Created: 07/05/2009
#
# Revisions:
# 07/05/2009    Sebastian Kelm    Moved from mescol.py into its own script, and generalised to all-atom data.
# 02/04/2012    Henry Wilman      Changed Residue   def __iter__(self), to make a CA optional

from array import array

import numpy

from prosci.common import join
from prosci.util.pdb import Pdb, residueLetter


class ResidueList(list):
  "A list of Residue objects."
  def __init__(self, pdb):
    "Takes a Pdb object as argument, or otherwise passes the argument to Pdb() first."
    if len(pdb) == 0:
      pass
    elif isinstance(pdb, ResidueList) or isinstance(pdb[0], Residue):
      for r in pdb:
        self.append(r)
    else:
      if not isinstance(pdb, Pdb):
        pdb = Pdb(pdb)
      for r in pdb.xresidues():
        self.append(Residue(r))
    try:
      self.code = pdb.code
    except:
      self.code = None
  
  def __repr__(self):
    return "ResidueList(%s)"%(list.__repr__(self))
  
  def __str__(self):
    return join("", self)
  
  def __getslice__(self, start=None, end=None):
    return ResidueList(list.__getslice__(self, start, end))
  
  def get_seq(self):
    "get_seq() : Returns a FASTA string representation of the sequence."
    seq = array('c')
    for r in self:
      seq.append(r.get_seq())
    return seq.tostring()
  
  def get_coords(self):
    "get a numpy array of coordinates"
    coords=[]
    for r in self:
      for a in r:
        coords.append(a.xyz)
    return numpy.array(coords)
  
  def deep_copy(self):
    rl = ResidueList([])
    for r in self:
      rl.append(Residue([a.copy() for a in r]))
    return rl
    


class Residue:
  "Class representing a single amino acid. Has fields for all main chain and CB atoms. Remaining atoms are in the list self.rest"
  def __init__(self, atoms):
    "Expects a list of Atom objects as input, or a Pdb object."

    self.type = None
    self.N = None
    self.CA = None
    self.C = None
    self.O = None
    self.CB = None
    self.rest = []
    
    for a in atoms:
      if a is None:
        continue
      self.type = a.res
      if a.atom == "N":
        self.N = a
      elif a.atom == "CA":
        self.CA = a
      elif a.atom == "C":
        self.C = a
      elif a.atom == "O":
        self.O = a
      elif a.atom == "CB":
        self.CB = a
      else:
        self.rest.append(a)
    
    if self.type in ("GLY", "ALA", "PRO"):
      del self.rest[:]
      if self.type == "GLY":
        self.CB = None

  
  def __iter__(self):
    if self.N is not None:
      yield self.N
    if self.CA is not None:
      yield self.CA
    if self.C is not None:
      yield self.C
    if self.O is not None:
      yield self.O

    if self.CB is not None:     
      yield self.CB             
    for a in self.rest:
      yield a
  
  
  def __str__(self):
    s=array('c')
    for a in self:
      s.extend(str(a))
      s.extend("\n")
    return s.tostring()
  
  
  def __repr__(self):
    return "Residue(%s)"%(repr([self.N, self.CA, self.C, self.O, self.CB]+self.rest))
  
  def get_seq(self):
    "Get the amino acid letter representing this Residue"
    return residueLetter(self.type)
  
  def set_type(self, newtype):
    """Set a new residues type.
    
    Deletes the side chain unless newtype == self.type. Keeps CB unless newtype == 'GLY'."""
    if newtype == self.type:
      return
    
    assert len(newtype) == 3
    assert newtype.isupper()
    
    if newtype == "GLY":
      self.CB = None
    del self.rest[:]
    
    self.type = newtype
    for a in self:
      a.res = newtype
