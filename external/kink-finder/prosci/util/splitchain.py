#!/usr/bin/python
# -*- coding: utf-8 -*-

import sys, os.path, re

from prosci.util.pdb import Pdb
from prosci.common import splitpath



re_default = re.compile("^(?:HEADER|TITLE|OBSLTE|COMPND|SOURCE|KEYWDS|EXPDTA|AUTHOR|REVDAT|SPRSDE|JRNL|REMARK|DBREF|SEQADV|FTNOTE|FORMUL|MODEL|ENDMDL|END   )")
re_atom    = re.compile("^(?:ATOM|HETATM|SHEET|TER)") # chainid at position 21
re_helix   = re.compile("^HELIX") # chainid at position 19
re_seqres  = re.compile("^SEQRES") # chainid at position 11
re_het     = re.compile("^HET    ") # chainid at position 12

def splitstructure(infile, prefix, chaincodes, extension=".pdb"):
  chaincodes = tuple(set(chaincodes))
  
  chainnum = range(len(chaincodes))
  
  if not (chaincodes and infile):
    return []
  
  outfiles = []
  outfnames = []
  for c in chaincodes:
    fname = prefix+c+extension
    outfiles.append(open(fname, 'w'))
    outfnames.append(fname)
  
  c = ' '
  m = None
  for line in infile:
    if re_default.match(line):
      for f in outfiles:
        f.write(line)
    else:
      m = re_atom.match(line)
      if m:
        c = line[21]
      else:
        m = re_helix.match(line)
        if m:
          c = line[19]
        else:
          m = re_seqres.match(line)
          if m:
            c = line[11]
          else:
            m = re_het.match(line)
            if m:
              c = line[12]
            else:
              continue
      for i in chainnum:
        if c == chaincodes[i]:
          outfiles[i].write(line)
          break
  
  for f in outfiles:
    f.close()
  
  return outfnames



def splitchains(files, options, doprint=False):
  for pdb_file in files:
    path, basename, ext = splitpath(pdb_file)
    
    pdb_code = basename
    
    a = Pdb(pdb_code, file(pdb_file))
    chains = a.get_chain_codes()
    
    if not chains or (len(chains)==1 and not chains[0]):
        #sys.stderr.write("No chain information found in PDB file '%s'. Not splitting it.\n" % (pdb_file))
        if ('f' in options) and ('a' in options):
            c = ''
            a_chain = a
            cgdb_id = ""
            if 'c' in options:
                cgdb_id = "CGDB{%s}" % (pdb_code.upper())
            text = ">%s\n%s%s\n%s\n" % (pdb_code+c, a_chain.get_structure_lign(), cgdb_id, a_chain.get_seq())
            f = open(basename+c+".ali", 'w')
            f.write(text)
            f.close()
            if os.path.isfile(basename+c+".ali"):
                if doprint:
                    print basename+c+".ali"
            else:
                sys.stderr.write("ERROR creating file: %s", basename+c+".ali")
    else:
        if 'p' in options:
                outfiles = splitstructure(file(pdb_file), basename, chains, ext)
                for f in outfiles:
                  if os.path.isfile(f):
                    if doprint:
                      print f
                  else:
                    sys.stderr.write("ERROR creating file: %s", f)
        
        if 'a' in options:
            for c in chains:
            #if 'p' in options:
                #os.system("cutchain %s %s > %s" % (c, pdb_file, basename+c+ext))
                #if os.path.isfile(basename+c+ext):
                    #if doprint:
                        #print basename+c+ext
                #else:
                    #sys.stderr.write("ERROR creating file: %s", basename+c+ext)
                
                a_chain = a.get_chain(c)
                
                #print "\n\n\n", str(a_chain), "\n\n\n"
                
                text = ">%s\n%s\n%s\n" % (pdb_code+c, a_chain.get_structure_lign(), a_chain.get_seq())
                f = open(basename+c+".ali", 'w')
                f.write(text)
                f.close()
                if os.path.isfile(basename+c+".ali"):
                    if doprint:
                        print basename+c+".ali"
                else:
                    sys.stderr.write("ERROR creating file: %s", basename+c+".ali")


if __name__ == "__main__":
    
    from prosci.shell import Params

    ################################
    # Command line option handling #
    ################################
    
    params = Params()
    
    
    if len(params.args) < 1 or ('p' not in params.opts and 'a' not in params.opts):
        print """
        Split a PDB or ATM file into chains.
        Output PDB structure files and/or ALI sequence files.
        
        USAGE:
            """+params.scriptname+""" OPTIONS pdb_file...
        
        SHORT OPTIONS:
            lower case letters enable an option  (e.g. -a   enables output in ALI format)
            UPPER CASE letters disable an option (e.g. -A  disables output in ALI format)
            
            OUTPUT FORMAT
            You must specify at least one of the following:
            
            -p     Write chains into PDB files (file extension of input files will be used)
            -a     Write chains into ALI files (.ali file extension will be used)
            
            -f     Force writing .ali files, even when only one chain is present.
            -c     Add CGDB entry name to structure lign, e.g. "CGDB{7AHL2G}"
        """
        sys.exit(1)
    
    
    ###############
    # File output #
    ###############
    
    splitchains(params.args, params.opts, doprint=True)
