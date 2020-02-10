# Benchmarks for Membrane Protein Energy Functions

This is a set of scientific benchmark tests to fit and evaluate membrane protein modeling energy functions. The test probe an energy function's ability to capture membrane protein orientation, stability, sequence, and structure. The methods are described in detail in the citation below. 

 - Alford RF & Gray JJ (2020) "Diverse scientific benchmarks reveal optimization imperatives for implicit membrane energy functions" _In Preparation_

## Manifest

 - `tests/` - Scientific benchmark scripts
 - `targets/` - Experimental and model datasets for testing
 - `analysis/` - Benchmark analysis scripts
 - `LICENSE` - MIT license for benchmark code
 - `config.txt` - Path and platform information for Rosetta

## Prerequisites

The test framework requires Python version 3.6 or 3.7. The tests also uses both the command-line and python interfaces to the Rosetta macromolecular modeling suite. Rosetta is available to academic users for free and to comercial users for a fee. 

To get Rosetta, obtain a license and download the package at <https://www.rosettacommons.org/software/license-and-download>. To compile the code, navigate to the `Rosetta/main/source/` directory and run the command below. 

```
./scons.py -j[XX] bin mode=release 
```

Here, "XX" is the number of processors you would like to compile with. For most users compiling on a laptop, the recommended number of processors is 1. If you are working on a larger workstation, we recommend scaling between 8-24 processors. More information can be found in the [Rosetta Build Documentation](https://www.rosettacommons.org/docs/wiki/build_documentation/Build-Documentation#setting-up-rosetta-3_basic-setup). 

To get PyRosetta, install miniconda first. On OSX or Linux, this is: 

```
wget https://repo.continuum.io/miniconda/Miniconda2-latest-MacOSX-x86_64.sh
bash Miniconda2-latest-MacOSX-x86_64.sh -b -p $HOME/miniconda
export PATH="$HOME/miniconda/bin:${PATH}""
```

You can then use conda to install PyRosetta. 

```
conda config --add channels "https://USERNAME:PASSWORD@conda.graylab.jhu.edu"
conda install pyrosetta
```

```
TODO: Include prerequisite of KinkFinder for helix kink calculations
Also include information about cluser ties for bigger tests (no laptop state)
```

Finally, the analysis code requires R version 3.0 or greater. You can download R here: https://www.r-project.org/. 

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

starting docs here. 


=== Everything below this line is in DRAFT and part of a TODO list ===

Drafting headers for python and R scripts

Add notes on cluster usage to this documentation as well as general reliance on supercomputers. 

Example for R
# @file: 	t01_tilt_angle.R
# @about: 	Map energies to orientations to find most stable configuration
# @author: 	Rebecca Alford <ralford3@jhu.edu>, Johns Hopkins University

Ideal interface

# Generate benchmark data
./generate_test_data.py --energy_fxn franklin2019 --which_tests all --restore_talaris False

	Option explanations
	--energy_fxn		Energy function to test, referred to by name of weights file in the Rosetta database. Must be present in both Rosetta & PyRosetta

	--which_tests		Include tests to run. Either all, or a comma separated list of one of the following tests given below: 

		TODO: add list, also refactor script to include labels 1-12




