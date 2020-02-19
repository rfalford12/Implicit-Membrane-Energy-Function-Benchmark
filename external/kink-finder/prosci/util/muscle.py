#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys
import os
#import tempfile
import subprocess
from array import array

from prosci.util.gaps import isGap, gappify


def reformat_seqs(seqs):
    if len(seqs) == 1:
      seqs = seqs[0]
    else:
      seqs = list(seqs)

    if not isinstance(seqs, str):
      for i in xrange(len(seqs)):
        if not isinstance(seqs[i], str):
          assert len(seqs[i]) == 2
          if not seqs[i][0].starswith(">"):
            seqs[i] = ">" + join("\n",seqs[i]).strip() + "\n"
          else:
            seqs[i] = join("\n",seqs[i]).strip() + "\n"
        elif not seqs[i].startswith(">"):
          seqs[i] = ">Seq%d\n%s\n" % (i, seqs[i].strip())
        elif not seqs[i].endswith("\n"):
          seqs[i] += "\n"
      seqs = "".join(seqs)
    return seqs


valid_residues = "ACDEFGHIKLMNPQRSTVWY-"
def replace_invalid_residues(txt):
    a=array('c')
    for line in txt.splitlines():
      if line.startswith(">"):
        a.extend(line)
      else:
        for c in line.strip():
          c=c.upper()
          if c in valid_residues:
            a.append(c)
          else:
            a.append('X')
      a.append("\n")
    return a.tostring()


def align_unstable(seqs):
    if not seqs:
      return ""
    
    p = subprocess.Popen("muscle -quiet", stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=sys.stderr, shell=True)
    text = p.communicate(seqs)[0]
    
    return text


def reorder_alignment(alignment_txt, seqs_txt):
    from prosci.util.fasta import Fasta
    
    alignment_txt = alignment_txt.rstrip()
    seqs_txt = seqs_txt.rstrip()
    assert (alignment_txt and seqs_txt) or not (alignment_txt or seqs_txt)
    if not alignment_txt:
      return ""
    
    alignment = Fasta("", alignment_txt) # unordered aligned sequences
    seqs = Fasta("", seqs_txt) # ordered unaligned sequences
    
    assert len(alignment) == len(seqs)
    
    # Index aligned sequences (speed-up for large alignments)
    alignment_entries = {}
    for e in alignment:
      alignment_entries[e.code] = e
    
    # Insert gaps into unaligned ordered sequences, while keeping original sequence
    for i in xrange(len(seqs)):
      seqs[i].seq = gappify(alignment_entries[seqs[i].code].seq, seqs[i].seq)
    
    # Return aligned ordered sequences
    return str(seqs)


def align(*seqs):
    seqs = reformat_seqs(seqs)
    if not seqs:
      return ""
    alignment = align_unstable(replace_invalid_residues(seqs))
    # Re-orders sequences AND makes sure input and output sequences have same letters (muscle replaces some things with X characters)
    #print alignment
    return reorder_alignment(alignment, seqs)
