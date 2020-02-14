# Benchmarks for Membrane Protein Energy Functions

This is a set of scientific benchmark tests for evaluating membrane protein modeling energy functions. The test probe an energy function's ability to capture membrane protein orientation, stability, sequence, and structure. The methods are described in detail in the citation below. 

 - Alford RF & Gray JJ (2020) "Diverse scientific benchmarks reveal optimization imperatives for implicit membrane energy functions" _In Preparation_

## Manifest

 - `tests/` - Scientific benchmark scripts
 - `targets/` - Experimental and model datasets for testing
 - `analysis/` - Benchmark analysis scripts
 - `LICENSE` - MIT license for benchmark code
 - `config.txt` - Path and platform information for Rosetta

## Prerequisites

#### System software

The test framework requires Python version 3.6+ and R version 3.6+. In addition, the data generation stage takes advantange of high-performance computing resources. the default setup is configured for a conda cluster. We also support SLURM clusters. For other setups, please email the author for assistance. 

#### Molecular modeling software

The tests also uses both the command-line and python interfaces to the Rosetta macromolecular modeling suite. Rosetta is available to academic users for free and to comercial users for a fee. 

To get Rosetta, obtain a license and download the package at <https://www.rosettacommons.org/software/license-and-download>. To compile the code, navigate to the `Rosetta/main/source/` directory and run the command below. 

```
./scons.py -jX bin mode=release 
```

Here, "X" is the number of processors to use during compilation. For compilation on a laptop, the recommended number of processors is 1. If you are working on a larger workstation or cluster, we recommend scaling between 8-24 processors. More information can be found in the [Rosetta Build Documentation](https://www.rosettacommons.org/docs/wiki/build_documentation/Build-Documentation#setting-up-rosetta-3_basic-setup). 

To get PyRosetta, install miniconda first. On OSX or Linux, this is: 

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:${PATH}""
```

You can then use conda to install PyRosetta, where USERNAME and PASSWORD are your Rosetta license credentials. 

```
conda config --add channels "https://USERNAME:PASSWORD@conda.graylab.jhu.edu"
conda install pyrosetta
```

We also use the KinkFinder package to compute helix kink angles. The KinkFinder package can be downloaded from [Charlotte Deane's lab website](http://opig.stats.ox.ac.uk/resources).

## Setup

To begin, modify the `config.txt` file to include information about your Rosetta installation, the location of this repository, and the platform. Each variable is described in the example below. 

```
benchmark = /path/to/benchmark  	# Path to this repo
rosettadir = /path/to/rosetta   	# Path to Rosetta bin
platform = linux 			# can be linux or mac
buildenv = release			# can be release or debug
compiler = gcc				# can be gcc or clang
```

## Documentation

The implicit membrane energy function benchmarks involves three steps: (1) data generation, (2) post-processing, and (3) analysis. We will walk through each step below. 

#### Step 1: Generate benchmark data

The first step performs all PyRosetta modeling calculations, any preparation steps for Rosetta applications, and submits calculations to the computing cluster. The generation script takes about 30min to run. Afterward, job completion requires roughly 1K CPU hours for all 12 tests. To run the generation step, use the command line below. 

	```
	./generate_test_data.py --energy_fxn franklin2019 --which_tests all
	Options:
	--energy_fxn		Energy function to test, referred to by name of weights 					file in the Rosetta database. Must be present in both 						Rosetta & PyRosetta
	--which_tests		Include tests to run. Either all, or a comma separated 						list of one of the following tests given below: 
	```

When the cluster jobs have completed, run the completion check script for step #1 to ensure that all of the data has been generated and there were no errors for individual jobs. To do so, use the command line below. 

	```
	./check_generate_step_complete.py --energy_fxn franklin2019 --which_tests all
	```

#### Step 2: Post-Process benchmark data

The sequence recovery and structure prediction benchmark tests require an imtermediate post-processing step before final data analysis. To run the post-processing step, run the command line below. 

	```
	./process_test_data.py --energy_fxn franklin2019 --which_tests all
	```

#### Step 4: 
