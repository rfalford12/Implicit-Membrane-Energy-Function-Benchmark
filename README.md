# Benchmarks for Membrane Protein Energy Functions

This is a set of scientific benchmark tests for evaluating membrane protein modeling energy functions. The test probe an energy function's ability to capture membrane protein orientation, stability, sequence, and structure. The methods are described in detail in the citation below. 

 - Alford RF & Gray JJ (2020) "Diverse scientific benchmarks reveal optimization imperatives for implicit membrane energy functions" _In Preparation_

## Manifest

 - `tests/` - Scientific benchmark scripts for generation, post-processing, and analysis
 - `targets/` - Experimental and model datasets for testing
 - `example/` - Example for analyzing output data from franklin2019
 - `LICENSE` - MIT license for benchmark code
 - `config.txt` - Path and platform information for Rosetta

## Prerequisites

#### System software

The test framework requires Python version 3.6+ and R version 3.6+. In addition, the data generation stage takes advantange of high-performance computing resources. the default setup is configured for a conda cluster. We also support SLURM clusters. For other setups, please email the author for assistance. 

#### Molecular modeling software

The tests use both the command-line and python interfaces to the Rosetta macromolecular modeling suite. Rosetta is available to academic users for free and to comercial users for a fee. 

To get Rosetta, obtain a license and download the package at <https://www.rosettacommons.org/software/license-and-download>. To compile the code, navigate to the `Rosetta/main/source/` directory and run the command below. 

```
./scons.py -jX bin mode=release 
```

Here, "X" is the number of processors to use during compilation. For compilation on a laptop, the recommended number of processors is 1. If you are working on a larger workstation or high performance computing cluster, we recommend scaling up to 8-24 processors. More information can be found in the [Rosetta Build Documentation](https://www.rosettacommons.org/docs/wiki/build_documentation/Build-Documentation#setting-up-rosetta-3_basic-setup). 

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

The first step performs all initial PyRosetta and Rosetta modeling calculations via a computing cluster. The generation script takes about 30min to run. Afterward, job completion requires roughly 1K CPU hours for all 12 tests. To run the generation step, use the command line below. 

	./generate_test_data.py --energy_fxn franklin2019 --which_tests all

The `--energy_fxn` flag sets the energy function to test, referred to by the name of the weights file in the Rosetta database (must be present in both Rosetta & PyRosetta). The `--which_tests` flag indicates which tests to be run. This can be `all` or a comma-separated list of the tests as given below. 

| Test                        | #  | Description 													   |
|-----------------------------|----|-------------------------------------------------------------------|
| tm-peptide-tilt-angle       | 1  | Tilt angle of single-span transmembrane peptides           	   |
| adsorbed-peptide-tilt-angle | 2  | Rotation angle of surface-adsorbed peptides    				   |
| protein-tilt-angle          | 3  | Tilt angle of multi-pass membrane proteins 					   |
| hydrophobic-length          | 4  | Hydrophobic thickness of multi-pass membrane proteins             |
| ddG-of-insertion            | 5  | Energetic cosf of transfering a peptide from water to bilayer     |
| ddG-of-pH-insertion         | 6  | Energetic cost of pH-dependent water-to-biolayer peptide transfer |
| ddG-of-mutation             | 7  | Energetic cost of single-point mutations in the membrane          |
| sequence-recovery           | 8  | Recovery of sequence features after full fixed-backbone redesign  |
| sc-distribution             | 9  | Depth-dependent distribution of side chains relative to native    |
| decoy-discrimination        | 10 | Native structure discrimination          						   |
| helix-kink                  | 11 | Helix kink angle prediction            						   |
| protein-protein-docking     | 12 | Membrane protein-protein docking           					   |

The outputs are then organized in a `data/` directory created by the script. The first subdirectory is the name of the energy function in use (e.g., franklin2019). Then, this folder contains 12 subdirectories for the data output from each test. 

#### Step 2: Post-Process benchmark data

The sequence recovery and structure prediction benchmark tests require an imtermediate post-processing step before final data analysis. To run the post-processing step, run the command line below. The flags are described in step #1. 

	./process_test_data.py --energy_fxn franklin2019 --which_tests all

#### Step 3: Analyze benchmark data 

The final step is to visualize and analyze the results of each benchmark tests. To do so, we provide a package of R scripts for generating the appropriate plots. While R can be run for the command line, we recommend downloading the R studio IDE from (https://rstudio.com/). You can run the example analysis script for `franklin2019` data from the `example/` directory through R studio or with the command line below. 

	Rscript analyze_f19_tests.R 

