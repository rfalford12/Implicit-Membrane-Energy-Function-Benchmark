# @file: plot_tilt_angle_and_depth_single_TM.R
# @author: Rebecca F. Alford (ralford3@jhu.edu)
# @brief: Compare experimental and predicted tilt angle and depth

library(nnet)
library(viridis)
library(ggplot2)
library(cowplot)

# For these - plot landscapes and both points
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/peptide-tilt-angle"
experimental.data <- read.table( paste( workdir, "experimental_tilt.dat", sep = "/"), header = T )
  
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

# Plot a grid of each energy landscape
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


# Read in orientation files and generate plot grids
df.tm.natives <- read.mapping(workdir, "A1_native_tm_ahelices")
p1 <- plot.orientation.map(df.tm.natives, "1a11", c(-45, -40, -35) )
p2 <- plot.orientation.map(df.tm.natives, "1mp6", c(-50, -40))
p3 <- plot.orientation.map(df.tm.natives, "1pje", c(-30, -25))
p4 <- plot.orientation.map(df.tm.natives, "2nr1", c(-40, -35, -30))
a1.grid <- plot_grid( p1, p2, p3, ncol = 3, nrow = 1, labels = c("a", "b", "c"), label_size = 9 )
#save_plot( paste( workdir, "SI_A1_native_tm_helices.pdf", sep = "/"), a1.grid, units = "in", base_width = 6, base_height = 2 )

df.surface.natives <- read.mapping(workdir, "A2_native_surface_ahelices")
q1 <- plot.orientation.map(df.surface.natives, "2mag", c(12, 16, 20) )
q2 <- plot.orientation.map(df.surface.natives, "1hu5", c(0, 2, 4, 6))
q3 <- plot.orientation.map(df.surface.natives, "1hu6", c(0, 2, 4))
q4 <- plot.orientation.map(df.surface.natives, "1hu7", c(-8, -6, -4))
q5 <- plot.orientation.map(df.surface.natives, "1f0d", c(32, 34, 36, 38))
q6 <- plot.orientation.map(df.surface.natives, "1f0e", c(-10, -8, -6))
q7 <- plot.orientation.map(df.surface.natives, "1f0g", c(35, 40, 45))
a2.grid <- plot_grid( q2, q3, q4, q5, q6, q7, labels = c("a", "b", "c", "d", "e", "f"), ncol = 3, nrow = 2 )
save_plot( paste( workdir, "SI_A2_native_surface_ahelices.pdf", sep = "/"), a2.grid, units = "in", base_width = 6, base_height = 4 )

df.tm.designs <- read.mapping(workdir, "A3_designed_tm_ahelices")
r1 <- plot.orientation.map(df.tm.natives, "WALP23")
r2 <- plot.orientation.map(df.tm.designs, "polyA-cappedW", c(-30, -26, -22) )
r3 <- plot.orientation.map(df.tm.designs, "polyA-cappedY")
a3.grid <- plot_grid(r1, r3, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 10 )
save_plot( paste( workdir, "SI_A3_designed_tm_ahelices.pdf", sep = "/"), a3.grid, units = "in", base_width = 4, base_height = 2 )

df.surface.designs <- read.mapping(workdir, "A4_designed_surface_ahelices")
s1 <- plot.orientation.map(df.surface.designs, "LK-peptide-n6", c(-45, -35))
save_plot( paste( workdir, "A4_designed_surface_ahelices.png", sep = "/"), s1, units = "in", base_width = 4, base_height = 4 )

# Make figure - four quadrants, 3/4 page figure
# 2nr1, polyA-cappedY, 2mag, lK-peptide n6
fig.grid <- plot_grid( p4,r2, q1, s1, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 10, scale = 0.97 )
save_plot( paste( workdir, "plots_peptide_orientation.pdf", sep = "/"), fig.grid, units = "in", base_width = 3.35, base_height = 3.5 )


