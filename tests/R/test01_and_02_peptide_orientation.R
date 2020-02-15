#' Analyze transmembrane and surface-adsorbed peptide orientation prediction
#' 
#' Computes analysis for peptide orientation calculations, including: 
#'    - single-span transmembrane peptides
#'    - helical surface-adsorbed peptides
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

library(nnet)
library(viridis)
library(ggplot2)
library(cowplot)


#' Compute the lowest energy orientation
#'
#' @param df Dataframe with mapping between peptide orientations and energies
#'
#' @return Dataframe with best tilt angle and depth position
#'
#' @examples
#' compute.best.orientation(df)
#' 
compute.best.orientation <- function( df ) { 

  ordered.df <- df[ order(df$total_score), ]
  head <- head(ordered.df,1)
  best.tilt <- head$angle
  best.depth <- head$zcoord

  # Make tilt "in-phase"
  if ( best.tilt > 270 ) {
    best.tilt = best.tilt - 270
  } else if ( best.tilt > 180 ) {
    best.tilt = best.tilt - 180
  } else if ( best.tilt > 90 ) {
    best.tilt <- best.tilt - 90
  }
  
  if ( best.tilt > 45 ) {
    best.tilt <- 90 - best.tilt
  }
  
  df <- data.frame( best.tilt = best.tilt, best.depth = best.depth )
  return(df)
}

#' Make a dataframe that maps energy landscapes with test information
#'
#' @param dir Working directory
#' @param test.name Name of benchmark test
#'
#' @return Dataframe with full energy landscapes mapped to target
#'
#' @examples
#' read.mapping( workdir, "A1" )
#' 
read.mapping <- function( dir, test.name ) {
  
  # Create a dataframe with the energy landscapes
  main.df <- data.frame()

  # List out all data files in the test category
  test.path = paste( dir, test.name, sep = "/" )
  data.files <- list.files(path = test.path, pattern = "*.dat")
  for(f in 1:length(data.files)) {
    
    # Grab the pdb identifier
    pdb.id <- strsplit(data.files[f], "_")[1]
    pdb.id <- unlist(pdb.id)[1]
    print(pdb.id)
    # Read data from landscape file
    file.name <- paste(test.path, data.files[f], sep = "/")
    df <- read.table( file.name, header = T )
    df$pdb.id <- pdb.id
    main.df <- rbind( main.df, df )
  }
  return(main.df)
}

#' Make a heatmap of orientations mapped to energies
#'
#' @param df Dataframe with orientations and energies
#' @param pdb.id PDB id for the arget of interest
#' @param score.breaks Range for colorbar
#'
#' @return Plot object
#'
#' @examples
#' plot.orientation.map(df, "1afo", c(-25, -5))
#' 
plot.orientation.map <- function(df, pdb.id, score.breaks) {

  # Get coordinates with the best orientation
  orientation.df <- compute.best.orientation(df[ which( df$pdb.id == pdb.id ), ])

  p <- ggplot() + 
    theme_bw() + 
    geom_raster( data = df[ which( df$pdb.id == pdb.id ), ], aes( x = angle, y = zcoord, fill = total_score ) ) + 
    geom_point( data = orientation.df, aes( x = best.tilt, y = best.depth ), color = "white" ) + 
    geom_point( data = orientation.df, aes( x = best.tilt, y = best.depth ), color = "white", shape = 1, size = 0.5 ) + 
    geom_point( data = experimental.data[ which( experimental.data$pdb == pdb.id ), ], aes( x = tilt, y = z ), color = "white" ) + 
    geom_point( data = experimental.data[ which( experimental.data$pdb == pdb.id ), ], aes( x = tilt, y = z ), color = "white", shape = 2, size = 2 ) + 
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(0, 90, 180, 270, 360) ) + 
    scale_y_continuous( "Depth (Å)", expand = c(0,0), limits = c(-45, 45) ) + 
    scale_fill_viridis( "Energy\n(REU)" ) + 
    theme( legend.position = "bottom",
           text = element_text( size = 7, color = "black"), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35, color = "black" ), 
           panel.border = element_rect(),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-2,-7))
  return(p)
}
