#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Common functions shared by all modules in prosci
#
# Author: Sebastian Kelm
# Created: 12/06/2009
#

import os
from array import array


def rxrange(start, end=None, step=1):
  assert step != 0, "rxrange() step must not be zero"
  if step > 0:
    if end==None:
      end=start
      start=0
    i=end-step
    while i >= start:
      yield i
      i -= step

def rrange(start, end=None, step=1):
  return list(rxrange(start, end, step))


def average(lst, ignoreNone=False):
  total=0.0
  n=0
  for x in lst:
    if not (ignoreNone and None==x):
      total += x
      n+=1
  return total/n

def AND(a, b):
  return a and b

def OR(a, b):
  return a or b

def NAND(a, b):
  return not (a and b)

def NOR(a, b):
  return not (a or b)

def XOR(a, b):
  return (a or b) and not (a and b)


def system(command):
  f = os.popen(command)
  text = f.read()
  f.close()
  return text


def system_list(command):
  f = os.popen(command)
  text = f.readlines()
  f.close()
  return text


# A generalised join function, that can join lists of non-string objects
# by Sebastian Kelm
#
def join(delim, lst):
    delim = str(delim)
    l = len(lst)
    if l < 1:
        return ""
    
    a = array('c')
    a.extend(str(lst[0]))
    i=1
    while i<l:
      a.extend(delim)
      a.extend(str(lst[i]))
      i+=1
    return a.tostring()


def stringKeepIndeces(source, indeces):
    #indeces.sort()
    a = array('c')
    for i in indeces:
      a.append(source[i])
    return a.tostring()


def stringRemoveIndeces(source, indeces):
    #indeces.sort()
    a = array('c')
    n=0
    for i in indeces:
      a.extend(source[n:i])
      n = i+1
    a.extend(source[n:len(source)])
    return a.tostring()


def listKeepIndeces(source, indeces):
    #indeces.sort()
    a = []
    for i in indeces:
      a.append(source[i])
    return a

def listKeepBoolean(source, keep):
    assert len(source) == len(keep)
    a = []
    for i in indeces:
      if keep[i]:
        a.append(source[i])
    return a


def listRemoveIndeces(source, indeces):
    #indeces.sort()
    a = []
    n=0
    for i in indeces:
      a.extend(source[n:i])
      n = i+1
    a.extend(source[n:len(source)])
    return a


def findAll(sequence, value):
    result=[]
    for i,v in enumerate(sequence):
      if v==value:
        result.append(i)
    return result


def findAllBoolean(sequence, value):
    result=[False] * len(sequence)
    for i,v in enumerate(sequence):
      if v==value:
        result[i] = True
    return result


def tokenize2(separators, seq):
  """tokenize2(separators, seq) : Transforms any type of sequence into a list of words and a list of starting indeces
  
  seq : the sequence (any sequence type, e.g. a list, tuple, string, ...) [will not be modified]
  separators : a sequence of values to be used as separators between words. adjoining separators will be merged.
  
  Returns: words, starting_indeces
           where
              len(starting_indeces) = len(words) + 1
              starting_indeces[-1] = len(seq)
  
  If seq is a string, words are also returned as strings. Otherwise every word is a list. Starting indeces are integers >= 0.
  """
  words = []
  starting_indeces = []
  if len(seq)<1:
    return words, [len(seq)]
  
  gapIsOpen=True # current gap size
  for i,v in enumerate(seq):
    if v not in separators:
      if gapIsOpen:
        words.append([])
        gapIsOpen=False
      words[-1].append(v)
    else:
      if not gapIsOpen:
        gaps.append(0)
        gapIsOpen=True
      gaps[-1] += 1
  
  if not gapIsOpen:
    gaps.append(0)
  
  assert len(gaps) == len(words) + 1
  
  if isinstance(seq, str):
    for i in xrange(len(words)):
      words[i] = "".join(words[i])
  
  starting_indeces.append(len(seq))
  return words, starting_indeces



def tokenize(separators, seq):
  """tokenize(separators, seq) : Transforms any type of sequence into a list of words and a list of gap lengths
  
  seq : the sequence (any sequence type, e.g. a list, tuple, string, ...) [will not be modified]
  separators : a sequence of values to be used as separators between words. adjoining separators will be merged.
  
  Returns: words, gaps
           where len(gaps) = len(words) + 1
  
  The first and last gaps are at the beginning and end of the sequence and may have length 0.
  
  If seq is a string, words are also returned as strings. Otherwise every word is a list. Gap lengths are integers >= 0.
  """
  words = []
  gaps  = [0]
  if len(seq)<1:
    return words, gaps
  
  gapIsOpen=True # current gap size
  for i,v in enumerate(seq):
    if v not in separators:
      if gapIsOpen:
        words.append([])
        gapIsOpen=False
      words[-1].append(v)
    else:
      if not gapIsOpen:
        gaps.append(0)
        gapIsOpen=True
      gaps[-1] += 1
  
  if not gapIsOpen:
    gaps.append(0)
  
  assert len(gaps) == len(words) + 1
  
  if isinstance(seq, str):
    for i in xrange(len(words)):
      words[i] = "".join(words[i])
  
  return words, gaps


def reduceTokens(words, gaps, wordcutoffs=0, minwordlength=1, wordValidator=None):
    """reduceTokens(words, gaps, wordcutoffs=0, minwordlength=1, wordValidator=None) : cuts off values from both sides of each token and deletes words that are then below a given length or that do not meet the condition given by the wordValidator(words, gaps, i) function. gap lengths are raised to reflect the changes. wordValidator should return True if the word is to be kept or False if the word should be deleted.
    
    Returns void. The input lists will be modified instead.
    """
    oldlength = sum(gaps)
    i=0
    while i < len(words):
      w = words[i]
      oldlength += len(w)
      if len(w) < (2*wordcutoffs)+minwordlength or (None != wordValidator and not wordValidator(words, gaps, i)):
        # Remove words that are too short
        gaps[i+1] += gaps[i] + len(w)
        del words[i]
        del gaps[i]
        continue
      else:
        # Modify words that are long enough
        words[i] = words[i][wordcutoffs:len(w)-wordcutoffs]
        gaps[i] += wordcutoffs
        gaps[i+1] += wordcutoffs
      i+=1
    
    newlength = sum(gaps)
    for w in words:
      newlength += len(w)
    
    assert newlength == oldlength


def deTokenize(separator, words, gaps):
    """deTokenize(separator, words, gaps) : reverse the work of the function tokenize(.). Returns a single sequence, where gaps are filled using separator.
    
    separator : a single value to be used to fill gaps
    words, gaps : lists of words and gap lengths, as returned by the function tokenize(.)
    """
    assert len(words) + 1 == len(gaps)
    seq=[]
    
    for i,w in enumerate(words):
      for j in xrange(gaps[i]):
        seq.append(separator)
      seq.extend(w)
    for j in xrange(gaps[-1]):
      seq.append(separator)
    
    if len(words) and isinstance(words[0], str):
      seq = "".join(seq)
    return seq


def binarySearch(lst, value, low=0, high=None):
    if high == None:
      high = len(lst)
    
    assert isinstance(low, int), "binarySearch : lower bound needs to be an integer. Given : "+repr(low)
    assert isinstance(high, int), "binarySearch : upper bound needs to be an integer. Given : "+repr(high)
    assert low <= high, "binarySearch : lower bound must be <= upper bound. Given: %d, %d" % (low, high)
    
    
    while low <= high:
        mid = (low + high) / 2
        if low == high:
          break
        if lst[mid] > value:
            high = mid
        elif lst[mid] < value:
            low  = mid + 1
        else:
            return mid
    
    return -1 - mid


def splitpath(somepath):
    # e.g.  /path/to/file.txt
    
    lastslash = somepath.rfind('/')
    
    path = somepath[:lastslash+1]
    filename = somepath[lastslash+1:]
    
    lastdot = filename.rfind('.')
    if lastdot<0:
        basename = filename
        ext = ""
    else:
        basename = filename[:lastdot]
        ext = filename[lastdot:]
    
    # e.g. "/path/to/", "file", ".txt"
    return (path, basename, ext)

def basename_noext(fname):
  return os.path.splitext(os.path.basename(fname))[0]


class ParsingError(RuntimeError):
    pass

class IllegalStateError(RuntimeError):
    pass
IllegalState = IllegalStateError

class DependencyError(RuntimeError):
    pass

class ArgumentError(ValueError):
    pass

class NotFoundError(IndexError):
    pass
