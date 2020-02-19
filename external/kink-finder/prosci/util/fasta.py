#!/usr/bin/python
# -*- coding: utf-8 -*-

from prosci.common import *
from prosci.util.gaps import gapchars
from prosci.util.seq import align


class Fasta(list):
  class Entry:
      def __init__(self, *args):
          if len(args) == 1:
            args = args[0]
          self.code, self.seq = args
          if self.code and self.code[0] == '>':
            self.code = self.code[1:]

      # Make Entry object seem like a list
      #
      def __iter__(self):
          yield self.code
          yield self.seq

      def __getitem__(self, i):
          if   0==i:
              return self.code
          elif 1==i:
              return self.seq
          raise IndexError("list index out of range")

      def __len__(self):
          return 2

      def __str__(self):
          return ">%s\n%s\n" % (self.code, self.seq)

      def __repr__(self):
          return "Fasta.Entry(%s,%s)" % (repr(self.code), repr(self.seq))

      def __cmp__(self, other):
          if isinstance(other, Fasta.Entry):
            return cmp(self.code, other.code)
          else:
            return cmp(self.code, str(other))

      def getCode(self):
          return self.code

      def get_seq(self):
          return self.seq

  def __init__(self, name, data=None):
    self.name = name
    if data:
      if isinstance(data, str):
        self.parseFastaFile(data.split("\n"))
      else:
        for e in data:
          self.append(Fasta.Entry(e))
    else:
      if isinstance(name, Fasta):
        self.name = name.name
        for e in name:
          self.append(Fasta.Entry(e))
      else:
        self.parseFastaFile(file(name))
  
  def __repr__(self):
    return "Fasta(%s)" % (list.__repr__(self))
  
  def __str__(self):
    return join("", self)
  
  def remove(self, code):
    n=-1
    for i,e in enumerate(self):
      if e.code == code:
        n=i
        break
    if n < 0:
      raise IndexError("Fasta entry not found: %s" %code)
    del self[n]
  
  def parseFastaFile(self, fileobj):
    for line in fileobj:
      line = line.strip()
      if not line:
        continue
      if line[0] == '>':
        if len(line) == 1:
          raise ParsingError("FASTA entry has zero-length title!")
        self.append(Fasta.Entry(line, ""))
      else:
        self[-1].seq += line
    if len(self) and not self[-1].seq:
      raise ParsingError("Last entry in FASTA file is empty!")

  def remove_gapped_columns(self, master=0):
    if isinstance(master, str):
      n = -1
      for i,e in enumerate(self):
        if e.code == master:
          n = i
          break
      if n < 0:
        print self
        raise KeyError("Entry not found in Fasta file: %s" %master)
      master = n
    
    residuesToDelete = []
    n = 0
    for i,c in enumerate(self[master].seq):
      if c in gapchars:
        residuesToDelete.append(i)

    for entry in self:
      assert len(entry.seq) == len(self[-1].seq), "Sequences not aligned"
      entry.seq = stringRemoveIndeces(entry.seq, residuesToDelete)
  
  def align(self):
    aligned = Fasta(self.name, align(str(self)))
    entries = {}
    for e in aligned:
      entries[e.code] = e
    
    assert len(aligned) == len(self)
    for i in xrange(len(self)):
      self[i] = entries[self[i].code]
    return self
