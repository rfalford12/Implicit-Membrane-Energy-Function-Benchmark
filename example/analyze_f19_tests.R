#' Analysis of results for scientific benchmarking of the franklin2019
#' all-atom implicit membrane energy function
#' 
#' Computes analysis for all scientific benchmark tests including: 
#'     - Test #1: Transmembrane peptide tilt angle
#'     - Test #2: Surface-adsorbed peptide rotation angle
#'     - Test #3: Multi-pass membrane protein orientation
#'     - Test #4: Hydrophobic Length
#'     - Test #5: Water-to-bilayer transfer energy for peptides
#'     - Test #6: pH-Dependent water-to-bilayer transfer energy for peptides
#'     - Test #7: Energetic cost of mutations
#'     - Test #8: Sequence Recovery
#'     - Test #9: Depth-dependent side chain distribution
#'     - Test #10: Decoy discrimination
#'     - Test #11: Helix Kinks
#'     - Test #12: Protein-protein docking
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

####################### Global Data #######################
benchmark.path <- "/home/ralford/research/Implicit-Membrane-Energy-Function-Benchmark"

###########################################################

source("../tests/R/test01_and_02_peptide_orientation.R")
source("../tests/R/test03_and_04_protein_orinetation.R")
source("../tests/R/test05_ddG_of_insertion.R")
source("../tests/R/test06_ddG_of_pH_insertion.R")
source("../tests/R/test07_ddG_of_mutation.R")
source("../tests/R/test08_sequence_recovery.R")
source("../tests/R/test09_side_chain_distribution.R")
source("../tests/R/test10_decoy_discrimination.R")
source("../tests/R/test11_helix_kinks.R")
source("../tests/R/test12_protein_docking.R")







