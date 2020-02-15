#' Analyze dG of peptide insertion results
#' 
#' Computes analysis for peptide insertion calculastions, including:
#'    - Correlation between predicted and experimentally measured value
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 
library(cowplot)
library(viridis)
library(ggrepel)

#' Compute the lowest energy orientation
#'
#' @param df Dataframe with mapping between peptide orientations and energies
#'
#' @return Dataframe with best tilt angle and depth position
#'
#' @examples
#' compute.best.orientation(df)
#'
compute.best.orientation <- function(df) {
  
  # Make a dataframe with only the membrane positions
  memb.df <- df[ which( abs(df$zcoord) <= 12 ), ]
  ordered.memb.df <- memb.df[ order(memb.df$total_score), ]
  head <- head(ordered.memb.df,1)
  best.memb.tilt <- head$angle
  best.memb.depth <- head$zcoord 
  best.memb.score <- head$total_score 

  # Make a dataframe with only the interfacial positions
  int.df <- df[ which( abs(df$zcoord) > 30  ), ]
  ordered.int.df <- int.df[ order(int.df$total_score), ]
  head <- head(ordered.int.df,1)
  best.int.tilt <- head$angle
  best.int.depth <- head$zcoord
  best.int.score <- head$total_score
  
  # Calculate the ddG of insertion
  ddG.of.insertion <- round( best.memb.score - best.int.score, 3 )
  
  # Create a dataframe that stores all of this information and return it
  out.df <- data.frame( best.memb.tilt = best.memb.tilt, 
                        best.memb.depth = best.memb.depth,
                        best.memb.score = best.memb.score, 
                        best.int.tilt = best.int.tilt, 
                        best.int.depth = best.int.depth, 
                        best.int.score = best.int.score, 
                        ddG.of.insertion = ddG.of.insertion )
  return(out.df)
}

#' Make a dataframe that maps energy landscapes with test information
#'
#' @param dir Working directory
#'
#' @return Dataframe with full energy landscapes mapped to target
#'
#' @examples
#' read.mapping( workdir )
#' 
read.mappings <- function( dir ) {
  
  main.df <- data.frame()
  data.dirs <- list.dirs(path = dir)
  for(d in 2:length(data.dirs)) {
    
    data.file <- list.files(path = data.dirs[d], pattern = "*.dat")
    pdb.id <- strsplit(data.file, "_")[1]
    file.name <- paste(data.dirs[d], data.file, sep = "/")

    df <- read.table( file.name, header = T )
    pdb.id <- unlist(pdb.id)[1]
    orientation.df <- compute.best.orientation(df)
    orientation.df$pdb.id <- pdb.id
    orientation.df$exp <- exp.data$ddG_insert[ which( exp.data$name == pdb.id ) ]
    orientation.df$error <- exp.data$error[ which( exp.data$name == pdb.id ) ]
    main.df <- rbind( main.df, orientation.df )
  }
  return(main.df)
}

#' Make correlation plot between predicted and experimentally measured values
#'
#' @param df Dataframe with predicted and experimentally measured values
#'
#' @return Plot Object
#'
#' @examples
#' make.dG.ins.correlation.plot(df)
#' 
make.dG.ins.correlation.plot <- function(df) {
  p <- ggplot( data = df, aes( x= exp, y = ddG.of.insertion, label = pdb.id ) ) + 
    theme_bw() + 
    background_grid() + 
    geom_vline( xintercept = 0, color = "gray50", linetype = "dashed", size = 0.35 ) + 
    geom_abline( slope = 2.396, intercept = -15.706, color = "gray50", size = 0.35 ) +
    geom_point( size = 0.75 ) + 
    geom_errorbarh( aes( xmax = exp+error, xmin = exp-error ) ) + 
    geom_text_repel( size = 2.5 ) + 
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-2.5, 2.5), expand = c(0,0)) + 
    scale_y_continuous( "Predicted (REU)", limits = c(-20, -10), expand = c(0,0) ) + 
    theme( text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35, color = "black" ), 
           panel.border = element_rect( color = "black" ))
  return(p)
}

