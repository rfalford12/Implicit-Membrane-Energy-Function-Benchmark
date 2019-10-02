#!/usr/bin/env python
###################################################################
#@file: 		predict_ideal_hydrophobic_length.py                                                                               
#@description: 	Compute the ideal hydrophobic length for a protein
#@args:			targetss
#@author: 		Rebecca F. Alford                   
#@email: 		rfalford12@gmail.com                                          
###################################################################

import random, sys, os
from pyrosetta import * 

def predict_hydrophobic_length( targets ): 
	"""
	A function to compute the energy minimized hydrophobic length
	Arguments: 
		targets = list of targets to evaluate
	"""

	# Steps of this script
	# 1. read all of the targets in as a long list of oriented PDB files
	# 2. add the membrane framework
	# 3. Backwards compute for t range 0-30 what 
	# 4. 
	# 5. 


