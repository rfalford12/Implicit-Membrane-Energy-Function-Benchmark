#' Analyze native structure discrimination (Test #10)
#' 
#' Computes analysis for native structure discrimination including: 
#'    - Score vs. RMSD plots for each target
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

library(ggplot2)
library(cowplot)
library(viridis)

#' Read refined decoy models into a data frame with score, target 
#' name, rms, and resolution
#' 
#' @param workdir Directory that includes scorefiles for all targets
#' 
#' @return Dataframe
#'
#' @example read.refined.models(workdir)
#' 
read.refined.models <- function( workdir ) {
  df <- data.frame()
  data.files <- list.files(path = workdir, pattern = "*.sc")
  for(f in 1:length(data.files)) {
    full.path <- paste( workdir, data.files[f], sep = "/" )
    sub.df <- read.table( full.path, header = T )
    sub.df$pdb <- unlist(strsplit( data.files[f], "_" )[1])[1]
    sub.df$resolution <- unlist(strsplit( data.files[f], "_" )[1])[2]
    sub.df <- subset( sub.df, select = c("total_score", "rms", "pdb", "resolution" ) )
    df <- rbind( df, sub.df )
  }
  return(df)
}

#' Add score and RMSD data for refined native models for each target
#' 
#' @param workdir Directory that includes refined_natives.txt
#' @param df Dataframe with score and rms values for each target
#' 
#' @return Dataframe
#'
#' @example add.refined.native.data(workdir, df)
#' 
add.refined.native.data <- function( workdir, df ) {
  refined.natives <- read.table( paste( workdir, "refined_natives.txt", sep = "/"), header = T )
  combined.df <- rbind( df, refined.natives )
  combined.df$resolution <- factor(combined.df$resolution, levels = c("alowres", "hires", "native"))
}

#' Make a score vs. rmsd plot for each target, colored by resolution
#' 
#' @param df Dataframe with score and rms values for each target
#' 
#' @return Plot object
#'
#' @example make.score.vs.rmsd.plot(df)
#' 
make.score.vs.rmsd.plot <- function(df) {
  p <- ggplot( data = combined.df, aes( x = rms, y = total_score, color = resolution ) ) +
    theme_bw() + 
    background_grid() + 
    geom_point( size = 0.35 ) +
    scale_x_continuous( "RMS (Å)", limits = c(0, 40) ) + 
    scale_y_continuous( "Score (REU)", limits = c(NA, -200 ) ) + 
    scale_color_manual( values = c("#377eb8", "#e41a1c", "#4daf4a") ) + 
    facet_wrap( ~ pdb, ncol = 5, nrow = 1, scales = "free_y") + 
    theme( legend.position = "none", 
           text = element_text( size = 8, color = "black" ), 
           axis.text = element_text( size = 8, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35 ), 
           panel.border = element_rect( color = "black" ))
  return(p)
    
}




