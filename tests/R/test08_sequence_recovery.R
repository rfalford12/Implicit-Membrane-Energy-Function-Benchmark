#' Analyze sequence recovery results (Test #8)
#' 
#' Computes analysis for sequence recovery, including: 
#'    - Recovery rates
#'    - KL Divergence
#'    - Amino acid composition
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

library(cowplot)
library(reshape2)
library(ggrepel)
library(viridisLite)
library(viridis)

#' Compute the fraction of amino acid types with recovery rates above random
#' or the general background probability of 0.05 (1 in 20)
#' 
#' @param df Dataframe with amino acid design and native distribution data
#'
#' @return Plot object
#'
#' @examples
#' frac.aa.recov.above.rand (d, "f19f)
#'
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

#' Make a plot to summarize the results of sequence recovery calculations 
#'
#' @param df Dataframe with amino acid design and native distribution data
#'
#' @return Plot object
#'
#' @examples
#' make.sequence.recovery.plot(df)
#'
make.sequence.recovery.plot <- function(df) {
  recovery.rates <- ggplot( data = df, aes( x = recovery, y = frac.recov.above.random*20, shape = solvation ) ) +
    background_grid() +
    geom_point(size = 1.5, color = "gray40") + 
    geom_text_repel( aes( label = efxn ), size = 2, force = 2 ) +
    scale_x_continuous( "Recovery (%)", limits = c(0, 0.5), expand = c(0,0) ) +
    scale_y_continuous( "# AAs with non-random Recovery", limits = c(0, 20), expand = c(0,0) ) +
    scale_shape_manual( values = c(15, 16, 1, 2, 17) ) + 
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
  return(recovery.rates)
}

#' Make a long-scale bar graph representing the KL divergence for each 
#' amino acid type from the native distribution for a given energy function's
#' design performance
#'
#' @param df Dataframe with amino acid design and native distribution data
#'
#' @return Plot object
#'
#' @examples
#' plot.aa.divergence(df)
#'
plot.aa.divergence <- function( df ) {
  p <- ggplot( df[ which( overall.df$variable == "overall" & df$efxn != "m19_all_DLPC" & df$efxn != "r15_all_DLPC" ), ], aes( x = residue, y = value, fill = category ) ) + 
    geom_hline( yintercept = 0, color = "black", size = 0.25 ) + 
    geom_bar( stat = "identity", color = "black", size = 0.15) + 
    scale_fill_brewer( palette = "Pastel1" ) + 
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
  return(p)  
}

#' Make a pie chart to represent tha distribution of amino acids in a set
#' of proteins, divided by amino acid chemical class
#'
#' @param df Dataframe with amino acid composition data
#'
#' @return Plot object
#'
#' @examples
#' make.aa.distribution.pie(df)
#'
make.aa.distribution.pie <- function(df) {
  p <- ggplot( df, aes( x = "", y = native, fill = Category ) ) + 
    geom_bar( stat = "identity" ) + 
    scale_fill_brewer( palette = "Pastel1" ) + 
    coord_polar("y", start = 0 ) + 
    theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )
  return(p)
} 
