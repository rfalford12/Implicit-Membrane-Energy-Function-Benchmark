#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# Hierarchical clustering algorithm
#
# Author: Sebastian Kelm
# Created: 25/04/2011
#


class TreeNode:
  def __init__(self, data=None, a=None, b=None):
    self.data = data
    self.a = a
    self.b = b
  
  def has_children(self):
    return self.a is not None or self.b is not None
  
  def count_leaves(self):
    if not self.has_children():
      return 1
    
    leaves = 0
    if self.a is not None:
      leaves += self.a.count_leaves()
    if self.b is not None:
      leaves += self.b.count_leaves()
    
    return leaves
  
  def get_first_leaf(self):
    x=self
    while x.a is not None:
      x = x.a
    return x.data
  
  def get_leaves(self):
    output_list = []
    self._get_leaves(output_list)
    return output_list
  
  def _get_leaves(self, output_list=None):
    if self.has_children():
      self.a._get_leaves(output_list)
      self.b._get_leaves(output_list)
    else:
      output_list.append(self.data)
  
  def __str__(self):
    if data is not None:
      return str(self.data)
    else:
      return "{"+str(self.a)+":"+str(self.b)+"}"
  
  def __repr__(self):
    return "TreeNode(%s, %s, %s)" % (repr(self.data), repr(self.a), repr(self.b))


def hierarchical_clustering(objects, simil_matrix, mode="single", cutoff=None):
    assert objects
    assert simil_matrix and len(simil_matrix) == len(simil_matrix[0])
    assert mode in ("single", "average"), "Illegal clustering mode: "+mode
    
    
    clusters = []
    for o in objects:
      clusters.append(TreeNode(o))
    
    
    # Keep merging the two most similar clusters, until there is only one left
    #
    while len(clusters) > 1:
    
      max_simil=0  # similarity of the two most similar clusters
      a=1  # matrix indeces of the most similar clusters
      b=0  # matrix indeces of the most similar clusters
      
      
      # Find the most similar pair of clusters
      #
      # NOTE: j is always smaller than i. Thus, b is always smaller than a
      #
      for i in xrange(1, len(clusters)):
        for j in xrange(i):
          if simil_matrix[i][j] > max_simil:
            max_simil = simil_matrix[i][j]
            a=i
            b=j
      
      
      if cutoff is not None and max_simil < cutoff:
        # Stop clustering, when we have clusters left that are less similar than a given cutoff
        break;
      
      
      # Merge the two clusters
      #
      newc = TreeNode(a=clusters[b], b=clusters[a])
      clusters[b] = newc
      del clusters[a]
      
      
      # Delete matrix rows & columns for the old clusters
      # Make a new matrix row & column for the merged cluster
      #
      if mode == "single":
          # Replace matrix row b with the "minimum distance" or "maximum similarity" values from rows a and b
          for i in xrange(len(simil_matrix)):
            simil_matrix[i][b] = simil_matrix[b][i] = Math.max(simil_matrix[b][i], simil_matrix[a][i])
      elif mode == "average":
          # Replace matrix row b with the "minimum distance" or "maximum similarity" values from rows a and b
          size_a = newc.a.count_leaves()
          size_b = newc.b.count_leaves()
          for i in xrange(len(simil_matrix)):
            simil_matrix[i][b] = simil_matrix[b][i] = (simil_matrix[b][i] * size_b + simil_matrix[a][i] * size_a) / (size_a + size_b)
      
      
      # Delete matrix row a, moving the following matrix rows up
      del simil_matrix[a]
      
      # Delete matrix column a, moving the following matrix columns left
      for i in xrange(len(simil_matrix)):
        del simil_matrix[i][a]
    
    
    assert clusters, "No clusters left... something went wrong"
    
    return clusters
