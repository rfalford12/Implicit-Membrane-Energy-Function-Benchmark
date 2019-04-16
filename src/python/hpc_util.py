#!/usr/bin/env python
###################################################################
#@file:         hpc_util.py                                                                                  
#@description:  Methods for submitting Rosetta jobs to SLURM and CONDOR clusters                                             
#@args:			jobname, executable, arguments, queue_no, high_mem                                                
#@author: 		Rebecca F. Alford                   
#@email: 		rfalford12@gmail.com                                          
###################################################################

import random, sys, os

def make_jobfile( casedir, case, executable, arguments ): 

    jobfile = casedir + "/" + case + "_energy_landscape.sh"
    with open( jobfile, 'a' ) as f: 
        f.write( "#!/bin/bash\n" )
        f.write( executable + " " + arguments + "\n" )
        f.close()
    os.system( "chmod +x " + jobfile )
    return jobfile


def submit_condor_job( path, jobname, executable, arguments, queue_no=1, high_mem=False ):

    # Given a test path and
    filename = path + "/" + jobname + ".condor"
    with open( filename, 'w' ) as f:
        f.write( "#!/bin/bash\n" )
        f.write( "# Automatically generated condor submission file for membrane benchmark\n" )
        f.write( "Universe = vanilla\n" )
        f.write( "output = " + path + "/" + jobname + ".out\n" )
        f.write( "error = " + path + "/" + jobname + ".err\n" )
        if ( high_mem == True ): 
            f.write( "request_memory = 15000\n")
        f.write( "notify_user = rfalford12@gmail.com\n" )
        f.write( "Executable = " + executable + "\n" )
        f.write( "Arguments = " + arguments + "\n" )
        f.write( "Queue " + str(queue_no) + "\n" )

    # Run the condor file
    os.system( "condor_submit " + filename )

def submit_stampede_job( path, jobname, jobfile, num_nodes=1 ): 
  
    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on Stampede2 with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH -J " + jobname + "\n" )
        f.write( "#SBATCH -p normal\n" )
        f.write( "#SBATCH -N " + str(num_nodes) + "\n" )
        f.write( "#SBATCH -n 16\n" )
        f.write( "#SBATCH -t 24:0:0\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH -o " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH -e " + path + "/" + jobname + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rfalford12@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "#SBATCH -A TG-MCB180056\n" )
        f.write( "\n" )

        # Specify required modules
        f.write( "module load intel/18.0.0\n" )
        f.write( "export MKL_MIC_ENABLE=1\n" )

        # Provide a description of the job
        f.write(  "echo Starting MPI job running " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpiexec bash " + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    #sbatch_command = "sbatch " + filename
    #os.system( sbatch_command )


def submit_marcc_job( path, jobname, jobfile, num_nodes=1 ): 

    # Create a new sbatch file named for the job type and test
    filename = path + "/" + jobname + ".sbatch" 
    with open( filename, 'w' ) as f: 

        # Write bash and comments
        f.write( "#!/bin/bash -l\n" )
        f.write ( "\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "# SLURM job script for membrane force field benchmarking applications\n" )
        f.write( "# Runs on MARCC with MPI applications\n" )
        f.write( "#----------------------------------------------------\n" )
        f.write( "\n" )

        # Write the job information
        f.write( "#SBATCH --job-name=" + jobname + "\n" )
        f.write( "#SBATCH --partition=parallel\n" )
        f.write( "#SBATCH --nodes=" + str(num_nodes) + "\n" )
        f.write( "#SBATCH --time=60:0:0\n" )
        f.write( "#SBATCH --mem=120GB\n" )

        # Write job specific output and reporting information
        f.write( "#SBATCH --output " + path + "/" + jobname + ".%j.out\n" )
        f.write( "#SBATCH --error " + path + "/" + jobname + ".%j.err\n" )
        f.write( "#SBATCH --mail-user=rfalford12@gmail.com\n" )
        f.write( "#SBATCH --mail-type=ALL\n" )
        f.write( "\n" )

        # Soecify required modules
        f.write( "module unload openmpi gcc\n" )
        f.write( "module load intel-mpi git\n" )
        f.write( "ml --gcc\n")
        f.write( "ml gcc/4.8.2\n")

        # Provide a description of the job
        f.write(  "echo Starting MPI job running " + jobfile + "\n" )

        # Run the job
        f.write( "time\n" )
        f.write( "mpirun ./" + jobfile + "\n" )
        f.write( "time\n" )

        f.close()

    # Run the slurm job file
    #sbatch_command = "sbatch " + filename
    #os.system( sbatch_command )
