
import os
import argparse
from Kink_Finder_Functions_Error import *

if __name__ == "__main__":
  
  
  # Make sure we can import stuff from this file's directory
  sys.path.append(os.path.abspath(os.path.dirname(sys.argv[0])))  
  
  #for use in finding the cylinder binary
  path = os.path.abspath(os.path.dirname(sys.argv[0]))
  
  parser = argparse.ArgumentParser()
  parser.add_argument('-f','--filename', action='store', 
                      default = 'all', help= 'Coordinate filename' +
                      ' e.g. 3EML.pdb',
                      metavar='')          
  parser.add_argument('-d', '--display', action= 'store_true', help = 'Display' +
                      ' summary on standard output')
  parser.add_argument('-o', '--output_dir', action= 'store', default = 'Output'+os.sep,
                     help = 'Directory for output',
                      metavar='')                    
  parser.add_argument('-p','--pymol', action= 'store_true',help = 'Write pymol files')
  parser.add_argument('-l','--helices', action= 'store', default = 'none', 
                      help = 'PDB number of first and last residues in each helix.'+ 
                      ' E.g "10-25 45-60" ',
                      metavar='') #does not work with more than one pdb_code
  
  args = parser.parse_args(sys.argv[1:])
  if len(sys.argv) == 1:
    parser.print_help() 
    exit()
  
  
  #=====================================================#
  # Arguments that can be changed from the command line #
  #=====================================================#
  directory = '.'+os.sep #where your pdb coordinate files are 
  #directory = '/data/cockatrice/wilman/data_sets/pdb_files_jan13/'
  
  user_helices = args.helices #only when running a single file
  pdb_extension = '.pdb' #what the file extensions of the pdb coordinate file is
  filename = args.filename #the pdb code of the file you wish to analyse. 
                # set this to 'all' to run every file in the directory
  
  display= args.display # display text output whilst running
  output_path = args.output_dir
  if args.pymol:
    pymol_file_directory = 'Pymol_files'+os.sep  #directory where the pymol files will be written to.
  else:
    pymol_file_directory = 'none'
  
  #output_path = 'Kink_Finder_results'+os.sep #where the result files will go
  if output_path[-1] != os.sep:
    output_path = output_path+os.sep #add a path seperator to the end of the output path, if it does not end with one
  
  #======================================================================#
  # More advanced arguments that cannot be changed from the command line #
  #======================================================================#  
  
  in_out = 'inside' #the kink residue is chosen to be on the inside of the helix 
  #in_out = 'outside' #the kink residue on the outside of the helix - uncomment if you wish to select kink residues on the outside of the kink
  
  break_angle = 100 # largest permissable angle in helix. Helices with angles larger than this are split in two, and the residue with the largest angle, and made into two seperate helices
  
  #======================================================================#
  # Other more involved arguments
  # These are for use with a .tem annotation file produced by Joy (K. Mizuguchi, C.M. Deane, T.L. Blundell, M.S. Johnson and J.P. Overington. (1998) JOY: protein sequence-structure representation and analysis. Bioinformatics 14:617-623), or DSSP (Dictionary of protein secondary structure: pattern recognition of hydrogen-bonded and geometrical features. Kabsch W, Sander C, Biopolymers. 1983 22 2577-2637). 
  # Changing the find_helices_from_tem variable in the kink_finder function within Kink_Finder_Functions.py lets kink finder find helices from a .tem file. 
  #======================================================================# 

  max_loop_length = 2 #if using a tem file, how many consecutive coil residues
                      # it is permissible to have within a helix
  soluble = True #ignoring the membrane, for use only with a tem file with membrane annotation - which can be obtained by using imembrane (Kelm S, Shi J, Deane CM. iMembrane: Homology-Based Membrane-Insertion of Proteins Bioinformatics. 2009). Setting to True, will, with a suitable annotation file, cause Kink Finder to only analyse helices that enter into the central tail region of the membrane.
  tem_extension = '.tem' #file extenstion of the annotation file. The annotation files are expected to have the same name as the structure files. E.g. 3EMLA.pdb would have annotation file 3EMLA.tem.
  
  kink_finder(directory, pdb_extension, tem_extension, filename, output_path, 
              soluble, display, break_angle, pymol_file_directory, 
              in_out, max_loop_length, user_helices, path)
