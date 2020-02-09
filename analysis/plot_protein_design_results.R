# @file: plot_protein_design_results.R
# @brief: Plot protein desgin results
# @author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(reshape2)
library(ggrepel)
library(scales)
library(dplyr)
library(viridisLite)
library(viridis)
library(expss)

# Current working directory
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/BioInspired-Membrane/figures/figure4-design-results-all"

# Files containing all of the data
total.seq.recov <- "total_sequence_recovery.dat"
total.kl.divergence <- "total_sequence_divergence.dat"
per.aa.seq.recov <- "per_amino_acid_recovery.dat"
per.aa.kl.divergence <- "per_amino_acid_divergence.dat"
per.class.seq.recov <- "per_amino_acid_class_recovery.dat"
per.class.kl.divergence <- "per_amino_acid_class_divergence.dat"

# For each energy function, calculate the fraction of per_AA recov that is above 0.05
per.aa.rates <- read.table( paste( workdir, per.aa.seq.recov, sep = "/"), header = T)
frac.aa.recov.above.rand <- function( per.aa.rates, efxn ) {
  
  overall.frac <- count_if( gt(0.05), per.aa.rates$overall[which( per.aa.rates$efxn == efxn)] )/20
  buried.frac <- count_if( gt(0.05), per.aa.rates$buried[which( per.aa.rates$efxn == efxn)] )/20
  surface.frac <- count_if( gt(0.05), per.aa.rates$surface[which( per.aa.rates$efxn == efxn)] )/20
  lipid_facing.frac <- count_if( gt(0.05), per.aa.rates$lipid_facing[which( per.aa.rates$efxn == efxn)] )/20
  interfacial.frac <- count_if( gt(0.05), per.aa.rates$interfacial[which( per.aa.rates$efxn == efxn)] )/20
  aqueous.frac <- count_if( gt(0.05), per.aa.rates$aqueous[which( per.aa.rates$efxn == efxn)] )/20
  
  # make a data frame with the info
  recovery.frac <- c( overall.frac, buried.frac, surface.frac, lipid_facing.frac, interfacial.frac, aqueous.frac )
  frac.labels <- c( "overall", "buried", "surface", "lipid_facing", "interfacial", "aqueous" ) 
  df <- data.frame( fraction = recovery.frac, solvation = frac.labels )
  df$efxn <- efxn
  return(df)
}

m07.recov.rates <- frac.aa.recov.above.rand( per.aa.rates, "m07_all_none" )
m12.recov.rates <- frac.aa.recov.above.rand( per.aa.rates, "m12_all_none" )
m19.recov.rates <- frac.aa.recov.above.rand( per.aa.rates, "m19_all_POPC" )
r15.recov.rates <- frac.aa.recov.above.rand( per.aa.rates, "r15_all_none" )

frac.recov.rates <- rbind( m07.recov.rates, m12.recov.rates, m19.recov.rates, r15.recov.rates )

# Load the overall sequence recovery
df.total.recov <- read.table( paste( workdir, total.seq.recov, sep = "/"), header = T )
trcov <- melt(df.total.recov, id.vars = c("efxn"), value.name = "recovery", variable.name = "solvation" )
tfrac <- frac.recov.rates[ order(frac.recov.rates$solvation),]

# # Some predefined target subsets
# all.targets <-  c( "m07_all_none", "m12_all_none", "r15_all_none", "m19_all_DLPC", "m19_all_POPC")

# load.design.performance.data <- function( recov.fn, div.fn, run.list ) {
#  
#   # Read the recovery and divergence data into melted dataframes
#   df.total.recov <- read.table( paste( workdir, recov.fn, sep = "/"), header = T )
#   trcov <- melt(df.total.recov, id.vars = c("efxn"), value.name = "recovery", variable.name = "solvation" )
#   df.total.div <- read.table( paste( workdir, div.fn, sep = "/"), header = T )
#   tdiv <- melt(df.total.div, id.vars = c("efxn"), value.name = "divergence", variable.name = "solvation" )
# 
#   # Merge into a single dataframe
#   total.stats <- trcov
#   total.stats$divergence <- tdiv$divergence
# 
#   # Use the "tag" to log the energy function version, species/taxonomy, and lipid types
#   efxn.versions <- c()
#   subset <- c()
#   lipid.types <- c()
#   for (efxn.type in total.stats$efxn) {
#     efxn.array <- strsplit(efxn.type, "_")
#     efxn.versions <- c( efxn.versions, efxn.array[[1]][1] )
#     subset <- c( subset, efxn.array[[1]][2] )
#     if ( length(efxn.array[[1]]) > 2 ) {
#       lipid.types <- c( lipid.types, efxn.array[[1]][3])
#     } else {
#       lipid.types <- c( lipid.types, "none")
#     }
#   } 
#   total.stats$efxn.version <- efxn.versions
#   total.stats$subset <- subset
#   total.stats$lipid.types <- lipid.types
#   
#   # Pull specific energy function versions for this data frame
#   total.stats.subset <- total.stats[ which( is.element( total.stats$efxn, run.list ) ), ]
#   return(total.stats.subset)
# }
# 
# plot.overall.design.performance <- function( df ) {
#   design.performance <- ggplot( data = df, aes( x = divergence, y = recovery, fill = lipid.types ) ) + 
#     background_grid() +
#     geom_point( aes( color = lipid.types ), size = 0.5 ) + 
#     geom_label_repel( aes( label = efxn.version ), force = 0.75, size = 2, label.padding = 0.125, label.size = 0.15, segment.size = 0.25 ) + 
#     geom_hline( yintercept = 0.05, color = "gray60", linetype = "dashed", size = 0.25 ) + 
#     geom_vline( xintercept = 0, color = "gray60", linetype = "dashed", size = 0.25 ) +
#     scale_y_continuous( "Sequence Recovery (%)", limits = c(0, NA) ) + 
#     scale_x_continuous( "KL Divergence", expand = c(0.1, 0.1) ) + 
#     scale_fill_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
#     scale_color_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
#     facet_wrap(  ~ solvation, ncol = 3, nrow = 2, scales = "free_x" ) + 
#     theme( legend.position = "none", 
#            text = element_text( size = 8 ), 
#            axis.text = element_text( size = 8 ), 
#            axis.line.x = element_line( size = 0.25 ), 
#            axis.line.y = element_line( size = 0.25 ), 
#            axis.ticks = element_line( size = 0.25 ), 
#            strip.background = element_blank() ) 
#   return(design.performance)
# }
# 
# # Plot the overall design performance for all five subsets
# all.targets.df <- load.design.performance.data( total.seq.recov, total.kl.divergence, all.targets )
# all.targets.overall.design <- plot.overall.design.performance( all.targets.df )
# save_plot( paste( workdir, "fig4_overall_design_performance.pdf", sep= "/"), all.targets.overall.design, units = "in", base_width = 3.42, base_height = 3 )
# 
# # Per AA design performance
# df.aa.recov <- read.table( paste( workdir, "per_amino_acid_recovery.dat", sep = "/"), header = T )
# aa.rcov <- melt(df.aa.recov, id.vars = c("efxn", "residue", "category"), value.name = "recovery", variable.name = "solvation" )
# #df.aa.div <- read.table( paste( workdir, "per_amino_acid_divergence.dat", sep = "/"), header = T )
# #aa.div <- melt(df.aa.div, id.vars = c("efxn", "residue", "category"), value.name = "kl.divergence", variable.name = "solvation" )
# per.aa.stats <- aa.rcov
# #per.aa.stats$kl.divergence <- aa.div$kl.divergence
# 
# design.all.list <-  c( "m07_all_none", "m12_all_none", "r15_all_none", "m19_all_POPC")
# per.aa.stats.all <- per.aa.stats[ which( is.element( per.aa.stats$efxn, design.all.list ) ), ]
# per.aa.design.perf.density <- ggplot( data = per.aa.stats.all, aes( x = recovery, fill = efxn, color = efxn ) ) +
#   background_grid() + 
#   geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) +
#   geom_density( alpha = 0.1 ) + 
#   scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5) ) +
#   scale_y_continuous( "Density", expand = c(0,0) ) +
#   scale_color_brewer(  palette = "Set1", direction = -1 ) + 
#   scale_fill_brewer( palette = "Set1", direction = -1 ) + 
#   facet_wrap( ~ solvation, scales = "free_y" ) +
#   theme_bw() + 
#   theme(  legend.position = "bottom",
#           axis.text = element_text( size = 8 ),
#           text = element_text( size = 8 ) )
# print(per.aa.design.perf.density)
# save_plot( paste( workdir, "SI_fig9_per_aa_design_distribution_lipid_facing.pdf", sep = "/"), per.aa.design.perf.density, units = "in", base_width = 6.5, base_height = 5 )
