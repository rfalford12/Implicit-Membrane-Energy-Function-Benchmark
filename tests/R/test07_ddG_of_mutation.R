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

substrRight <- function(x, n){
  substr(x, nchar(x)-n+1, nchar(x))
}

# Read in the ddG breakdown
breakdown.df <- read.table( paste( workdir, "C3_OmpLA_aro_ddGs", "breakdown2.sc", sep = "/"), header = T )
breakdown.trp.df <- breakdown.df[ which( str_sub(breakdown.df$mutation, -1,-1) == "W"), ]
breakdown.tyr.df <- breakdown.df[ which( str_sub(breakdown.df$mutation, -1,-1) == "Y"), ]

trp.df <- melt( breakdown.trp.df, ids = c("mutation") ) 
tyr.df <- melt( breakdown.tyr.df, ids = c("mutation") )

p <- ggplot( data = trp.df, aes( x = variable, y = value ) ) + 
  theme_bw() + 
  background_grid() + 
  geom_hline( yintercept = 0, size = 0.35 ) + 
  geom_bar( position = "dodge", stat = "identity" ) + 
  scale_x_discrete( "", expand = c(0,0) ) + 
  scale_y_continuous( "Score (REU)" ) + 
  facet_wrap( ~ mutation, ncol = 4, scales = "free_x" ) + 
  theme( text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.text.x = element_text( angle = 90, hjust = 1 ),
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35, color = "black" ), 
         panel.border = element_rect( color = "black" ))
save_plot( paste( workdir, "trp_breakdown.pdf", sep = "/"), p, units = "in", base_width = 6.8, base_height = 7 )

q <- ggplot( data = tyr.df, aes( x = variable, y = value ) ) + 
  theme_bw() + 
  background_grid() + 
  geom_hline( yintercept = 0, size = 0.35 ) + 
  geom_bar( position = "dodge", stat = "identity" ) + 
  scale_x_discrete( "", expand = c(0,0) ) + 
  scale_y_continuous( "Score (REU)" ) + 
  facet_wrap( ~ mutation, ncol = 4, scales = "free_x" ) + 
  theme( text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.text.x = element_text( angle = 90, hjust = 1 ),
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35, color = "black" ), 
         panel.border = element_rect( color = "black" ))
save_plot( paste( workdir, "tyr_breakdown.pdf", sep = "/"), q, units = "in", base_width = 6.8, base_height = 7 )


df3 <- read.table( paste( workdir, "C3_OmpLA_aro_ddGs", "breakdown3.sc", sep = "/"), header = T )

plot.ddG.of.mutation.by.depth.v2 <- function( df, test.name ) {
  
  # Grab list of mutations with experimental data
  p <- ggplot( data = df[ which( df$mutation != "P" & df$mutation != "F" & df$mutation != "W" ), ], aes( x = exp, y = dG_memb, fill = abs(depth), label = Pos) ) + 
    theme_bw() + 
    background_grid() + 
    geom_hline( yintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_vline( xintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_abline( size = 0.35 ) + 
    geom_label( size = 2.5, label.padding = unit(0.15, "lines") ) + 
    scale_fill_viridis( "Depth (Å)", breaks = c(0, 4, 8, 12), begin = 0.2 ) + 
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-5.5, 2.5), expand = c(0,0) )  +
    scale_y_continuous( "Predicted (REU)", limits = c(-4, 5), expand = c(0,0) ) + 
    theme( legend.position = "bottom",
           text = element_text( size = 8 ), 
           axis.text = element_text( size = 8, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35 ), 
           panel.border = element_rect( color = "black" ),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-7,-7))
  
  q <- ggplot( data = df[ which( df$mutation != "P" & df$mutation != "F" & df$mutation != "Y" ), ], aes( x = exp, y = dG_memb, fill = abs(depth), label = Pos) ) + 
    theme_bw() + 
    background_grid() + 
    geom_hline( yintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_vline( xintercept = 0, color = "gray40", size = 0.35 ) + 
    geom_abline( size = 0.35 ) + 
    geom_label( size = 2.5, label.padding = unit(0.15, "lines") ) + 
    scale_fill_viridis( "Depth (Å)", breaks = c(0, 4, 8, 12), begin = 0.2 ) + 
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-5.5, 2.5), expand = c(0,0) )  +
    scale_y_continuous( "Predicted (REU)", limits = c(-5, 5), expand = c(0,0) ) + 
    theme( legend.position = "bottom",
           text = element_text( size = 8 ), 
           axis.text = element_text( size = 8, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35 ), 
           panel.border = element_rect( color = "black" ),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-7,-7))
  
  y.and.w <- plot_grid( p,q, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 10)
  
  return(y.and.w)
}
