#!/usr/bin/env python
###################################################################
#@file:         read_config.py                                                                                  
#@description:  Read configuration file, output benchmark & rosetta path                                             
#@args:			--energy_fxn (wts) --targets (list of targets)                                                          
#@author: 		Rebecca F. Alford                   
#@email: 		rfalford12@gmail.com                                          
###################################################################

import sys, os

configpath = "../../config.txt"

class BenchmarkConfig(): 
	def __init__( benchmark_path, rosetta_path, platform, buildenv, compiler ): 
		self.benchmark_path = benchmark_path
		self.rosetta_path = rosetta_path
		self.platform = platform
		self.buildenv = buildenv
		self.compiler = compiler

def read_config(): 

	with open( configpath, 'rt' ) as f: 
		contents = f.readlines()
		contents = [ x.strip() for x in contents ]
		contents = [ x.split(" ") for x in contents ]

	config = BenchmarkConfig( contents[0][2], contents[1][2], contents[2][2], contents[3][2], contents[4][2] )
	return config
