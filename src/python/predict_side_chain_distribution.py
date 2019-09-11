#!/usr/bin/env python
###################################################################
#@file:         predict_side_chain_distributions.py                                
#@description:  Calculate the depth-dependent side chain 
#					distributions from designed PDBS                                                       
#@author:     	Rebecca F. Alford                   
#@email:    	rfalford12@gmail.com                                          
###################################################################

import random, os, sys
import hpc_util, read_config
from string import Template
from progress.bar import Bar
from pyrosetta import * 

init()

class AminoAcidInfo: 
	def __init__(self): 
		self.coords = []

def compute_side_chain_distribution( config, native_pdbs, designed_pdbs ): 
	"""
		A function for calculating the depth dependent amino acid distributions

		Arguments: 
			config = Path to benchmark and Rosetta executables
			native_pdbs = List of native PDBs
			designed_pdbs = List of designed PDBs
	"""

	# Read pdblist and write to src file for natives
	native_poses = read_pdbs( native_pdbs )
	get_side_chain_distribution( native_poses, "native" )

	# Read pdblist and write to src file for designed pdbs
	designed_poses = read_pdbs( designed_pdbs )
	get_side_chain_distribution( designed_poses, "design" )

def read_pdbs( pdbfile_list ): 

	with open( pdbfile_list, 'rt' ) as f: 
		contents = f.readlines()
		contents = [ x.strip() for x in contents ]
	f.close()

	pdbs = []
	for pdbfile in contents: 
		if ( os.path.isfile( pdbfile ) ): 
			pdb = pose_from_pdb( pdbfile )
			pdbs.append( pdb )
		else: 
			sys.exit( "PDB File " + pdbfile + " not found!" )
	
	return pdb

def get_side_chain_distribution( pdblist, src ): 

	print("Generating distributions for source " + src )

	# Initialize an AA info object for each amino acid
	ala = AminoAcidInfo()
	cys = AminoAcidInfo()
	asp = AminoAcidInfo()
	glu = AminoAcidInfo()
	phe = AminoAcidInfo()
	gly = AminoAcidInfo()
	his = AminoAcidInfo()
	ile = AminoAcidInfo()
	lys = AminoAcidInfo()
	leu = AminoAcidInfo()
	met = AminoAcidInfo()
	asn = AminoAcidInfo()
	pro = AminoAcidInfo()
	gln = AminoAcidInfo()
	arg = AminoAcidInfo()
	ser = AminoAcidInfo()
	thr = AminoAcidInfo()
	val = AminoAcidInfo()
	trp = AminoAcidInfo()
	tyr = AminoAcidInfo()

	# Make a progress bar
	bar = Bar('Reading amino acids in PDB files', max=len(pdblist))

	# Loop through each PDB and encode case statement
	for pdb in pdblist: 
		for i in range(1, pdb.total_residue()+1): 

			if ( pdb.residue(i).name1() == "A" ): 
				ala.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "C" ): 
				cys.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "D" ): 
				asp.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "E" ): 
				glu.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "F" ): 
				phe.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "G" ): 
				gly.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "H" ): 
				his.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "I" ): 
				ile.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "K" ): 
				lys.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "L" ): 
				leu.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "M" ):
				met.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "N" ): 
				asn.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "P" ): 
				pro.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "Q" ): 
				gln.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "R" ): 
				arg.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "S" ): 
				ser.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "T" ): 
				thr.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "V" ): 
				val.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "W" ): 
				trp.coords.append( pdb.residue(i).xyz("CA").z() )
			elif ( pdb.residue(i).name1() == "Y" ): 
				tyr.coords.append( pdb.residue(i).xyz("CA").z() )

		bar.next()
	bar.finish()

	# Write each of the distributions to an output data file... 
	filename = src + "_side_chain_distribution.txt"
	with open( filename, 'wt' ) as f: 

		# write file header
		f.write( "AA zcoord src" )

		# Write for each coordinate list
		aa = "A"
		print("Writing coordinates for alanine...")
		for coord in ala.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "C"
		print("Writing coordinates for cystine...")
		for coord in cys.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "D"
		print("Writing coordinates for aspartate...")
		for coord in asp.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "E"
		print("Writing coordinates for glutamate...")
		for coord in glu.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "F"
		print("Writing coordinates for phenylaalanine...")
		for coord in phe.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "G"
		print("Writing coordinates for glycine...")
		for coord in gly.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "H"
		print("Writing coordinates for histidine...")
		for coord in his.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "I"
		print("Writing coordinates for isoleucine...")
		for coord in ile.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "K"
		print("Writing coordinates for lysine...")
		for coord in lys.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "L"
		print("Writing coordinates for leucine...")
		for coord in leu.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "M"
		print("Writing coordinates for methionine...")
		for coord in met.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "N"
		print("Writing coordinates for asparagine...")
		for coord in asn.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "P"
		print("Writing coordinates for proline...")
		for coord in pro.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "Q"
		print("Writing coordinates for glutamate...")
		for coord in gln.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "R"
		print("Writing coordinates for arginine...")
		for coord in arg.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "S"
		print("Writing coordinates for serine...")
		for coord in ser.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "T"
		print("Writing coordinates for threonine...")
		for coord in thr.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "V"
		print("Writing coordinates for valine...")
		for coord in val.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "W"
		print("Writing coordinates for tryptophan...")
		for coord in trp.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		aa = "Y"
		print("Writing coordinates for tyrosine...")
		for coord in tyr.coords: 
			f.write( f'{aa} {coord} {src}\n' )

		f.close()
