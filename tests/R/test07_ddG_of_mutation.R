#' Analyze ddG of single-point mutation results (Test #7)
#' 
#' Computes analysis for ddG of mutation calculastions, including:
#'    - Correlation between experimentally measured and predicted values
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 
library(cowplot)
library(viridis)
library(reshape2)
library(ggplot2)
library(stringr)

#' Make a plot of predicted vs. experimentally measured ddG of mutation and
#' label the colors by amino acid category
#'
#' @param workdir Working directory for the stability measurement results
#' @param test.name Name of test dataset (e.g. C2_PagP_canonical_ddGs)
#'
#' @return Plot Object
#'
#' @examples
#' plot.ddG.of.mutation.by.AA( workdir, "C1_OmpLA_canonical_ddGs")
#'
plot.ddG.of.mutation.by.AA <- function( dir, test.name ) {
  
  # Grab list of mutations with experimental data
  test.path = paste( dir, test.name, sep = "/" )
  data.file <- list.files( path = test.path, pattern = "*.dat")
  file.path = paste( test.path, data.file, sep = "/")
  df <- read.table( file.path, header = T )
  
  p <- ggplot( data = df[ which( df$Mut != "P" ), ], aes( x = experimental_ddG, y = predicted_ddG, fill = class, label = Mut) ) + 
    background_grid() + 
    geom_hline( yintercept = 0, color = "gray40" ) + 
    geom_vline( xintercept = 0, color = "gray40" ) + 
    geom_abline() + 
    geom_label() + 
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-5.5, 5.5), expand = c(0,0) ) + 
    scale_y_continuous( "Predicted (REU)", limits = c(-5.5, 5.5), expand = c(0,0) ) 
    
  return(p)
}

#' Make a plot of predicted vs. experimentally measured ddG of mutation and
#' label the colors by amino acid depth relative to the membrane center
#'
#' @param workdir Working directory for the stability measurement results
#' @param test.name Name of test dataset (e.g. C2_PagP_canonical_ddGs)
#'
#' @return Plot Object
#'
#' @examples
#' plot.ddG.of.mutation.by.depth( workdir, "C1_OmpLA_canonical_ddGs")
#'
plot.ddG.of.mutation.by.depth <- function( dir, test.name ) {
  
  # Grab list of mutations with experimental data
  test.path = paste( dir, test.name, sep = "/" )
  data.file <- list.files( path = test.path, pattern = "*.dat")
  file.path = paste( test.path, data.file, sep = "/")
  df <- read.table( file.path, header = T )
  
  p <- ggplot( data = df[ which( df$Mut != "P" & df$Mut != "F" & df$Mut != "W" ), ], aes( x = experimental_ddG, y = predicted_ddG, fill = abs(depth), label = Pos) ) + 
    theme_bw() + 
    background_grid() + 
    geom_hline( yintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_vline( xintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_abline( size = 0.35 ) + 
    geom_label( size = 2.5, label.padding = unit(0.15, "lines") ) + 
    scale_fill_viridis( "Depth (Å)", breaks = c(0, 4, 8, 12), begin = 0.2 ) + 
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-5.5, 2.5), expand = c(0,0) )  +
    scale_y_continuous( "Predicted (REU)", limits = c(-5, 12), expand = c(0,0) ) + 
    theme( legend.position = "bottom",
           text = element_text( size = 8 ), 
           axis.text = element_text( size = 8, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35 ), 
           panel.border = element_rect( color = "black" ),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-7,-7))
  
  q <- ggplot( data = df[ which( df$Mut != "P" & df$Mut != "F" & df$Mut != "Y" ), ], aes( x = experimental_ddG, y = predicted_ddG, fill = abs(depth), label = Pos) ) + 
    theme_bw() + 
    background_grid() + 
    geom_hline( yintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_vline( xintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_abline( size = 0.35 ) + 
    geom_label( size = 2.5, label.padding = unit(0.15, "lines") ) + 
    scale_fill_viridis( "Depth (Å)", breaks = c(0, 4, 8, 12), begin = 0.2 ) + 
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-5.5, 2.5), expand = c(0,0) )  +
    scale_y_continuous( "Predicted (REU)", limits = c(-5, 12), expand = c(0,0) ) + 
    theme( legend.position = "bottom",
           text = element_text( size = 8 ), 
           axis.text = element_text( size = 8, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35 ), 
           panel.border = element_rect( color = "black" ),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-7,-7))
  
  y.and.w <- plot_grid( p,q, ncol = 2, nrow = 1, labels = c("b", "c"), label_size = 10)
  
  return(y.and.w)
}

