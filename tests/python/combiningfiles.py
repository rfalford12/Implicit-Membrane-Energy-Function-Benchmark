#!/usr/bin/env python
""" Master script for combininig all files into one file
and generating a energy landscape as a function of depth and tilt 
angle  minimized  over rotation (/azimuthal) angle. 
This module is used for tests 1,2,3,5 and 6. 

Authors: 
	Rituparna Samanta <rsamanta@utexas.edu> 

Example: 
	$ python3 combiningfiles.py --which_tests "tm-peptide-tilt-angle" --energy_fxn "franklin2019"

Arguments: 
	- energy_fxn: Weights file for energy function of interest
	- which_tests: Run all tests, or specify by comma-separated list

Requirements: 
	- Rosetta release 246 or greater
	- PyRosetta4 for Python 3.6 or 3.7
"""

import sys, os
import read_config
import glob
import math
import fileinput
import pandas as pd
import numpy as np
from optparse import OptionParser, IndentedHelpFormatter
_script_path_ = os.path.dirname( os.path.realpath(__file__) )

def unique(list1):
    unique_list = []
    for x in list1:
        if x not in unique_list:
            unique_list.append(x)
    return(unique_list)

def combiningallfiles(peptide_list):
     
    for element_num in range(len(peptide_list)):
        
        start_index = ['-60.000000','-40.000000','-20.000000', '0.000000', '20.000000','40.000000']
        end_index = ['-40.000000','-20.000000', '0.000000', '20.000000', '40.000000','60.000000']#[-3, -2, -1, 0, 0, 0, 1, 2]
          
        peptide_name = peptide_list[element_num]

        if( peptide_name!= "LK_peptide_n6"):
            peptide_name = peptide_name.split('_')[0]

        labelin1 = peptide_list[element_num] + '/' + peptide_name + '_combined.dat'
        labelout1 = peptide_list[element_num] + '/' + peptide_name + '_min.dat'
        print(labelin1)
        print(labelout1)
        
          
         
        label0= peptide_list[element_num]+'/'+peptide_name+'_franklin2019_'+str(end_index[0])+'_'+str(start_index[0])+'_landscape.dat'
        df0 = pd.read_csv(label0)    
        label1= peptide_list[element_num]+'/'+peptide_name+'_franklin2019_'+str(end_index[1])+'_'+str(start_index[1])+'_landscape.dat'
        df1 = pd.read_csv(label1)
        label2= peptide_list[element_num]+'/'+peptide_name+'_franklin2019_'+str(end_index[2])+'_'+str(start_index[2])+'_landscape.dat'
        df2 = pd.read_csv(label2)
        label3= peptide_list[element_num]+'/'+peptide_name+'_franklin2019_'+str(end_index[3])+'_'+str(start_index[3])+'_landscape.dat'
        df3 = pd.read_csv(label3)
        label4= peptide_list[element_num]+'/'+peptide_name+'_franklin2019_'+str(end_index[4])+'_'+str(start_index[4])+'_landscape.dat'
        df4 = pd.read_csv(label4)
        label5= peptide_list[element_num]+'/'+peptide_name+'_franklin2019_'+str(end_index[5])+'_'+str(start_index[5])+'_landscape.dat'
        df5 = pd.read_csv(label5)
        

        df = pd.concat([df0,df1,df2,df3,df4,df5])
        df.to_csv(labelin1, index=False)
        df_read = pd.read_csv(labelin1, delimiter=" ")
        X = unique(df_read['zcoord']) 
        Y = unique(df_read['angle'])
        
        mindata = []

        for i in range(len(X)):
            
            for j in range(len(Y)):
                
                newarr = df_read[df_read['angle']==Y[j]]
                secondarr = newarr[newarr['zcoord']==X[i]]
                arr = np.array(secondarr['total_score']) 
                minpos = np.argmin(arr)
                mindata.append(secondarr.iloc[minpos,0:23])
        mindatadf = pd.DataFrame(mindata)        
        mindatadf.to_csv(labelout1, sep=" ", index=False)        
        
        del df0, df1, df2
        del df3, df4, df5
        del df, df_read
        del X, Y 
    
   # fclose(fileID);


def main( args ):       
    
    all_sub_tests = ["ddG-of-insertion","ddG-of-pH-insertion","tm-peptide-tilt-angle","adsorbed-peptide-tilt-angle","protein-tilt-angle"]

    # Read options from the command line
    parser = OptionParser( usage="usage %prog --energy_fxn franklin2019 --which_tests all" )
    parser.set_description(main.__doc__)
	
    parser.add_option( '--energy_fxn', '-e', action="store", help="Name of energy function weights file", )
    parser.add_option( '--which_tests', '-w', action="store", help="Specify tests run (comma separated list)", )

    (options, args) = parser.parse_args(args=args[1:])
    global Options
    Options = options

    # Check that required options have been provided
    if ( not Options.energy_fxn or not Options.which_tests ): 
       print("Missing required options --which_tests or --energy_fxn" )
       sys.exit()
    
    # Check test categories
    test_names = []

    # Read path configuration file
    config = read_config.read_config()
    
    if (  Options.which_tests == "all" ): 
       test_names = all_sub_tests
    else: 
       test_names = Options.which_tests.split(",")
       # check that all names are valid
       for name in test_names: 
           if name not in all_sub_tests: 
                  sys.exit( "No such test " + name + ". Exiting! or not among the tests 1-3 and 5-6" )
            
           else:
                  
                 
                  datadir = config.benchmark_path + "data/" + Options.energy_fxn + "/" + name + "/"
                  print(datadir)
                  if ( not os.path.isdir(datadir) ): 
                       sys.exit( "No such test data available; test needs to finish before combining files" )
                  os.chdir( datadir )

                  if( name == "tm-peptide-tilt-angle" ):
                       peptide_list = ['1a11','1mp6','WALP23','2nr1','1pje','polyA-cappedW','polyA-cappedY']
                  elif( name == "adsorbed-peptide-tilt-angle" ):
                       peptide_list = ['1f0d','1f0e','1f0g','1hu5','1hu6','1hu7','2mag','LK_peptide_n6']
                  elif( name == "protein-tilt-angle" ):
                       peptide_list = ['1fep','1gzm','1m0l','1nqe','1okc','1p4t','1qfg','1qj8','1qjp','1r3j','1yce','2cfp','2j8c','2qom','2x9k','3aeh','3dzm','3pox','3syb','3wbn','3wfd','4afk','4fqe','4hyj','4m48','4n6h','4rl8','4uc2','4x5n']
                  elif( name == "ddG-of-insertion"):
                       peptide_list = ['GL5','GL6','GL7','GL8','GWL6']
                  elif( name == "ddG-of-pH-insertion"):
                       peptide_list = ['pHLIP-v1_4','pHLIP-v2_8','pHLIP-v2_4','pHLIP-v2_8','pHLIP-v3_4','pHLIP-v3_8','pHLIP-v4_4','pHLIP-v4_8','pHLIP-v5_4','pHLIP-v5_8','pHLIP-v6_4','pHLIP-v6_8','pHLIP-v7_4','pHLIP-v7_8','pHLIP-v8_4','pHLIP-v8_8','pHLIP-v9_4','pHLIP-v9_8','pHLIP-v10_4','pHLIP-v10_8','pHLIP-v11_4','pHLIP-v11_8','pHLIP-v12_4','pHLIP-v12_8','pHLIP-v13_4','pHLIP-v13_8','pHLIP-v14_4','pHLIP-v14_8','pHLIP-v15_4','pHLIP-v15_8','pHLIP-v16_4','pHLIP-v16_8']
                
                    
                  combiningallfiles(peptide_list)

if __name__ == "__main__" : main(sys.argv)            
