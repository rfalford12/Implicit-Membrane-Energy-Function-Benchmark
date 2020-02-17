# @file: plot_protein_design_results.R
# @brief: Plot protein desgin results
# @author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(reshape2)
library(ggrepel)
library(viridisLite)
library(viridis)


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


frac.recov.rates <- rbind( m07.recov.rates, m12.recov.rates, m19.recov.rates, r15.recov.rates )


load.design.performance.data <- function( recov.fn, div.fn, run.list ) {

  # Read the recovery and divergence data into melted dataframes
  df.total.recov <- read.table( paste( workdir, recov.fn, sep = "/"), header = T )
  trcov <- melt(df.total.recov, id.vars = c("efxn"), value.name = "recovery", variable.name = "solvation" )
  df.total.div <- read.table( paste( workdir, div.fn, sep = "/"), header = T )
  tdiv <- melt(df.total.div, id.vars = c("efxn"), value.name = "divergence", variable.name = "solvation" )

  # Merge into a single dataframe
  total.stats <- trcov
  total.stats$divergence <- tdiv$divergence

  # Use the "tag" to log the energy function version, species/taxonomy, and lipid types
  efxn.versions <- c()
  subset <- c()
  lipid.types <- c()
  for (efxn.type in total.stats$efxn) {
    efxn.array <- strsplit(efxn.type, "_")
    efxn.versions <- c( efxn.versions, efxn.array[[1]][1] )
    subset <- c( subset, efxn.array[[1]][2] )
    if ( length(efxn.array[[1]]) > 2 ) {
      lipid.types <- c( lipid.types, efxn.array[[1]][3])
    } else {
      lipid.types <- c( lipid.types, "none")
    }
  }
  total.stats$efxn.version <- efxn.versions
  total.stats$subset <- subset
  total.stats$lipid.types <- lipid.types

  # Pull specific energy function versions for this data frame
  total.stats.subset <- total.stats[ which( is.element( total.stats$efxn, run.list ) ), ]
  return(total.stats.subset)
}

plot.overall.design.performance <- function( df ) {
  design.performance <- ggplot( data = df, aes( x = divergence, y = recovery, fill = lipid.types ) ) +
    background_grid() +
    geom_point( aes( color = lipid.types ), size = 0.5 ) +
    geom_label_repel( aes( label = efxn.version ), force = 0.75, size = 2, label.padding = 0.125, label.size = 0.15, segment.size = 0.25 ) +
    geom_hline( yintercept = 0.05, color = "gray60", linetype = "dashed", size = 0.25 ) +
    geom_vline( xintercept = 0, color = "gray60", linetype = "dashed", size = 0.25 ) +
    scale_y_continuous( "Sequence Recovery (%)", limits = c(0, NA) ) +
    scale_x_continuous( "KL Divergence", expand = c(0.1, 0.1) ) +
    scale_fill_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
    scale_color_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
    facet_wrap(  ~ solvation, ncol = 3, nrow = 2, scales = "free_x" ) +
    theme( legend.position = "none",
           text = element_text( size = 8 ),
           axis.text = element_text( size = 8 ),
           axis.line.x = element_line( size = 0.25 ),
           axis.line.y = element_line( size = 0.25 ),
           axis.ticks = element_line( size = 0.25 ),
           strip.background = element_blank() )
  return(design.performance)
}

per.aa.design.perf.density <- ggplot( data = per.aa.stats.all, aes( x = recovery, fill = efxn, color = efxn ) ) +
  background_grid() +
  geom_vline( xintercept = 0.05, color = "gray60", linetype = "dashed" ) +
  geom_density( alpha = 0.1 ) +
  scale_x_continuous( "Sequence Recovery (%)", limits = c(0, 0.5) ) +
  scale_y_continuous( "Density", expand = c(0,0) ) +
  scale_color_brewer(  palette = "Set1", direction = -1 ) +
  scale_fill_brewer( palette = "Set1", direction = -1 ) +
  facet_wrap( ~ solvation, scales = "free_y" ) +
  theme_bw() +
  theme(  legend.position = "bottom",
          axis.text = element_text( size = 8 ),
          text = element_text( size = 8 ) )

# Current working directory
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/BioInspired-Membrane/figures/figure4-design-results-all"
df <- read.table( paste( workdir, "design_seqrecov_performance.txt", sep = "/"), header = T )
df$solvation <- factor( df$solvation, c("overall", "buried", "surface", "lipid_facing", "interfacial", "aqueous"))
df$which_panel <- factor( df$which_panel, c("FALSE", "burial", "memb"))
recovery.rates <- ggplot( data = df, aes( x = recovery, y = frac.recov.above.random*20, shape = solvation ) ) +
  background_grid() +
  geom_point(size = 1.5, color = "gray40") + 
  geom_text_repel( aes( label = efxn ), size = 2, force = 2 ) +
  scale_x_continuous( "Recovery (%)", limits = c(0, 0.5), expand = c(0,0) ) +
  scale_y_continuous( "# AAs with non-random Recovery", limits = c(0, 20), expand = c(0,0) ) +
  scale_shape_manual( values = c(15, 16, 1, 2, 17) ) + 
  #scale_fill_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
  #scale_color_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
  facet_wrap(  ~ which_panel, ncol = 3, nrow = 1, scales = "free" ) +
  theme_bw() + 
  theme( legend.position = "none",
         text = element_text( size = 8 ),
         axis.text = element_text( size = 8 ),
         axis.line.x = element_line( size = 0.25 ),
         axis.line.y = element_line( size = 0.25 ),
         axis.ticks = element_line( size = 0.25 ),
         strip.background = element_blank(), 
         panel.grid.major = element_line( size = 0.15 ), 
         panel.grid.minor = element_line( size = 0.15 ), 
         panel.border = element_rect( size = 0.25 ) ) + 
  coord_flip()

p <- ggplot( overall.df[ which( overall.df$variable == "overall" & overall.df$efxn != "m19_all_DLPC" & overall.df$efxn != "r15_all_DLPC" ), ], aes( x = residue, y = value, fill = category ) ) + 
  geom_hline( yintercept = 0, color = "black", size = 0.25 ) + 
  geom_bar( stat = "identity", color = "black", size = 0.15) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  #scale_fill_manual( values = c( "#984ea3", "#e41a1c", "#ffff33", "#377eb8", "#ff7f00") ) + 
  scale_y_continuous( limits = c(-5, 1), breaks = c(-5, -4, -3, -2, -1, 0, 1), "Divergence" ) + 
  scale_x_discrete( "" ) + 
  facet_wrap( ~ efxn, scales = "free" ) + 
  theme_bw() + 
  theme( strip.background = element_blank(), 
         legend.position = "none", 
         text = element_text( size = 8 ), 
         axis.text = element_text( size = 7 ), 
         axis.line.x = element_line( size = 0.25 ), 
         axis.line.y = element_line( size = 0.25 ), 
         axis.ticks = element_line( size = 0.25 ), 
         panel.grid.major = element_line( size = 0.15 ), 
         panel.grid.minor = element_line( size = 0.15 ), 
         panel.border = element_rect( size = 0.25 ) )
print(p)

q <- ggplot( overall.df[ which( overall.df$variable == "lipid_facing" & overall.df$efxn != "m19_all_DLPC" ), ], aes( x = residue, y = value, fill = category ) ) + 
  geom_hline( yintercept = 0, color = "black" ) + 
  geom_bar( stat = "identity", color = "black", size = 0.15) + 
  scale_fill_manual( values = c( "#984ea3", "#e41a1c", "#ffff33", "#377eb8", "#ff7f00") ) + 
  scale_y_continuous( limits = c(-16, 5), "Divergence" ) + 
  scale_x_discrete( "" ) + 
  facet_wrap( ~ efxn, scales = "free" ) + 
  theme_bw() + 
  theme( strip.background = element_blank(), 
         legend.position = "none", 
         text = element_text( size = 8 ), 
         axis.text = element_text( size = 7 ), 
         axis.line.x = element_line( size = 0.25 ), 
         axis.line.y = element_line( size = 0.25 ), 
         axis.ticks = element_line( size = 0.25 ), 
         panel.grid.major = element_line( size = 0.25 ), 
         panel.grid.minor = element_line( size = 0.25 ) )
print(q)

# Plot the native distribution
native.pie <- ggplot( df2, aes( x = "", y = native, fill = Category ) ) + 
  geom_bar( stat = "identity" ) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )

m07.pie <- ggplot( df2, aes( x = "", y = m07, fill = Category ) ) + 
  geom_bar( stat = "identity") + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )

m12.pie <- ggplot( df2, aes( x = "", y = m12, fill = Category ) ) + 
  geom_bar( stat = "identity" ) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )

m19.pie <- ggplot( df2, aes( x = "", y = m19, fill = Category ) ) + 
  geom_bar( stat = "identity" ) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank() )

