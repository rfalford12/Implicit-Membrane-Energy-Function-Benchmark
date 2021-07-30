#' Analyze multi-pass membrane protein orientation and hydrophobic thickness
#' 
#' Computes analysis for membrane protein orientation calculations, including: 
#'   - Cummulative distribution of residuals for tilt, depth, and thickness
#'   - Comparison of OPM reference and predicted values
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu> and updated by RituparnaSamanta <rsamant2[at]jhu.edu>
#' 

library(ggplot2)
library(viridis)
library(cowplot)
library(nnet)
library(ggrepel)


#' Compute the lowest energy orientation
#'
#' @param df Dataframe with mapping between peptide orientations and energies
#'
#' @return Dataframe with best tilt angle and depth position
#'
#' @examples
#' orientation(df)
#'
compute.best.orientation <- function( df ) {
  #min.index <- which.min( df$fa_water_to_bilayer )
  min.index <- which.min( df$total_score )
  best.tilt <- df$angle[ min.index ]
  best.depth <- df$zcoord[ min.index ]
  best.score <- df$total_score[ min.index ]

  df <- data.frame( best.tilt = best.tilt, best.depth = best.depth, best.score = best.score )
  return(df)
}


#' Make a dataframe that maps energy landscapes with test information
#'
#' @param dir Working directory, experiment file name
#'
#' @return Dataframe with full energy landscapes mapped to target
#'
#' @examples
#' read.mapping( workdir )
#'
#' 

read.mappings <- function( dir, experimentalfilename ) {

  main.df <- data.frame()
  opm.orientation <- read.table( paste( dir, experimentalfilename, sep = "/"), header = T )
 
  data.dirs <- list.dirs(path = dir)
  
  for(d in 2:length(data.dirs)) {

    data.file <- list.files(path = data.dirs[d], pattern = "*.dat")
  
    pdb.id <- unlist(strsplit(data.file[1], "_"))[1]
    
    print(pdb.id)
    
    sub.df <- data.frame()
    
   
    file.name <- paste(data.dirs[d], paste(pdb.id,"combined.dat", sep = "_"), sep = "/")
    sub.df <-read.table( file.name, header = T)
    #colnames(sub.df) = c("zcoord", "angle", "azimuthal", "total_score", "fa_water_to_bilayer")
    
    #extracting the minimum over azimuthal angle data frame
    sub_min.df <- orientation.frame(data.dirs[d], sub.df, pdb.id, 0)
    #extracting the dataframe at azimuthal angle = 0
    sub_prev.df <- orientation.frame(data.dirs[d], sub.df, pdb.id,  1)
    #sub_min.df <- sub_prev.df
    
    p1 <- make.energy.plots(sub_min.df)
    p2 <- make.energy.plots(sub_prev.df)
    plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
    a1.grid <-plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
    #a1.grid <-plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a"), label_size = 9)
    save_plot( paste( data.dirs[d], paste(pdb.id,sep="_","full.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )
    
    orientation.df <- compute.best.orientation(sub_min.df)
    orientation_prev.df <- compute.best.orientation(sub_prev.df)
    
    orientation.df$tilt_prev <- orientation_prev.df$best.tilt
    orientation.df$depth_prev <- orientation_prev.df$best.depth
    orientation.df$score_prev <- orientation_prev.df$best.score
    orientation.df$pdb.id <- pdb.id
    orientation.df$alpha_or_beta <- opm.orientation$alpha_or_beta[ which( opm.orientation$pdb == pdb.id ) ]
     
    orientation.df$opm.tilt <- opm.orientation$angle[ which( opm.orientation$pdb == pdb.id ) ]
    orientation.df$opm.depth <- opm.orientation$depth[ which( opm.orientation$pdb == pdb.id ) ]
    orientation.df$opm.delcenter <- opm.orientation$delcenter[ which( opm.orientation$pdb == pdb.id ) ]

    main.df <- rbind( main.df, orientation.df )
    
  }
  
  a <- main.df$best.tilt[ which(main.df$best.tilt > 90 ) ]
  b <- (180-a)
  main.df$best.tilt[ which(main.df$best.tilt >90 ) ] <- b
  
  a <- main.df$best.tilt[ which(main.df$best.tilt < -90 ) ]
  b <- (-180-a)
  main.df$best.tilt[ which(main.df$best.tilt < -90 ) ] <- b
  
  a <- main.df$tilt_prev[ which(main.df$tilt_prev > 90 ) ]
  b <- (180-a)
  main.df$tilt_prev[ which(main.df$tilt_prev >90 ) ] <- b
  
  a <- main.df$tilt_prev[ which(main.df$tilt_prev < -90 ) ]
  b <- (-180-a)
  main.df$tilt_prev[ which(main.df$tilt_prev < -90 ) ] <- b
  
  
  write.table( main.df, file = paste( dir, "processed_data.txt", sep = "/") )
  
  
  return(main.df)
}

read.mappings.hl <- function( dir, experimentalfilename) {
  
  main.df <- data.frame()
  opm.orientation <- read.table( paste( workdir, experimentalfilename, sep = "/"), header = T )
  
  data.dirs2 <- list.dirs(path = dir)
  
  for(d in 2:length(data.dirs2)) {
    
    data.file <- list.files(path = data.dirs2[d], pattern = "*thickness.dat")
    
    pdb.id <- unlist(strsplit(data.file[1], "_"))[1]
    
    print(pdb.id)
    
    sub.df <- data.frame()
    
    # for( index in 2:length(data.file) ){ 
    #  file.name <- paste(data.dirs[d], data.file[index], sep = "/")
    # print(file.name)
    #  df <- read.table( file.name, header = T )
    #  sub.df <- rbind(sub.df,df)
    # }
    file.name <- paste(data.dirs2[d], data.file, sep = "/")
    sub.df <-read.table( file.name, header = T)
    
    min.index <- which.min( sub.df$score )
    best.thk <- sub.df$thickness[ min.index ]
    
    orientation.frame <- data.frame( opm.thk = opm.orientation$length[ which( opm.orientation$pdb == pdb.id ) ], best.thk = best.thk, pdb = pdb.id, alpha_or_beta = opm.orientation$alpha_or_beta[ which( opm.orientation$pdb == pdb.id ) ] )
    main.df <- rbind( main.df, orientation.frame )
  }
  write.table( main.df, file = paste( dir, "processed_data_thk.txt", sep = "/") )
  
  
  return(main.df)
}

#' Compute and output residuals relative to the reference value
#'
#' @param dir Dataframe of best tilt angle and depth values
#'
#' @return None - outputs a file called "protein_position.dat" with the data
#'
#' @examples
#' compute.and.output.residuals(df)
#' 
compute.and.output.residuals <- function( df ) {
  a <- df$angle[ which(df$angle > 90 ) ]
  b <- (180-a)
  df$angle[ which(df$angle >90 ) ] <- b
  
  a <- df$angle[ which(df$angle < -90 ) ]
  b <- (-180-a)
  df$angle[ which(df$angle < -90 ) ] <- b
  
  
  #Residue wrt to y=x curve.
  # have taken absolute value of the depths, since the system is symmetric about the 
  #center of the membrane and 180 degrees. a protein at x degrees and depth = y is same as at 180-x degrees and depth = -y. 
  #Rosetta calculated center of the protein for rotation whereas OPM of the TM section. 
 
   df$depth.res <- ( abs(df$opm.depth) - abs(df$best.depth + df$opm.delcenter) )
  
  df$depth.res_prev <- ( abs(df$opm.depth) - abs(df$depth_prev + df$opm.delcenter) )
 
  df$angle.res <- ( df$opm.tilt - df$best.tilt )

  df$angle.res_prev <- ( df$opm.tilt - df$tilt_prev )

  write.table( df, paste( workdir, "protein_position.dat", sep = "/"))
  return(df)
  
}

#' Make plot to compare predicted and reference tilt angle values
#'
#' @param df Dataframe with predicted and OPM reference values
#'
#' @return Plot object
#'
#' @examples
#' make.tilt.angle.comparison.plot(df)
#' 
make.tilt.angle.comparison.plot <- function(df, flag) {
  
  for.ribbon <- seq(-90, 90, by=1)
 
  
  p <- ggplot() +
    theme_bw() + 
    background_grid() + 
    geom_abline( size = 0.35 ) +
    geom_ribbon( aes(ymin=for.ribbon-10, ymax=for.ribbon+10, x=for.ribbon), alpha = 0.2) +
    if(flag==0){
      geom_point( data = df, aes( x = opm.tilt, y = best.tilt, color = alpha_or_beta ), size = 0.75 )  
    }else{
      geom_point( data = df, aes( x = opm.tilt, y = tilt_prev, color = alpha_or_beta ), size = 0.75 ) 
    } 
   q<- p + coord_cartesian(ylim = c(0, 90), xlim =c(0, 90) ) +
    scale_x_continuous( "OPM (degrees ˚)" ) +
    scale_y_continuous( "Predicted (˚)" ) + 
    scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
    theme( legend.position = "top", text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           panel.border = element_rect( color = "black"), 
           axis.ticks = element_line( size = 0.35 ) )
    
  return(q)
}

#' Make plot to compare predicted and reference depth values
#'
#' @param df Dataframe with predicted and OPM reference values
#'
#' @return Plot object
#'
#' @examples
#' make.depth.comparison.plot(df)
#' 
make.depth.comparison.plot <- function(df, flag) {
  for.ribbon <- seq(-50, 50, by=1)
  q <- ggplot() + 
    theme_bw() + 
    background_grid() + 
    geom_abline( size = 0.35 ) + 
    geom_ribbon( aes(ymin=for.ribbon-2.5, ymax=for.ribbon+2.5, x=for.ribbon), alpha = 0.2) +
    #if(flag==0){
    #  geom_point( data = df, aes( x = abs(opm.depth), y = abs(best.depth + opm.delcenter), color = alpha_or_beta ), size = 0.75  ) 
    #}else{
    #  geom_point( data = df, aes( x = abs(opm.depth), y = abs(depth_prev + opm.delcenter), color = alpha_or_beta ), size = 0.75 )  
    #}
  if(flag==0){
    geom_point( data = df, aes( x = abs(opm.depth), y = abs(best.depth), color = alpha_or_beta ), size = 0.75  ) 
  }else{
    geom_point( data = df, aes( x = abs(opm.depth), y = abs(depth_prev), color = alpha_or_beta ), size = 0.75 )  
  }
   
  p<- q +  coord_cartesian(ylim = c(0, 20), xlim =c(0, 20) ) +
    scale_x_continuous( "OPM (Å)" ) +
    scale_y_continuous( "Predicted (Å)" ) +
    scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
    theme( legend.position = "top",
           text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           panel.border = element_rect( color = "black"), 
           axis.ticks = element_line( size = 0.35 ) )
  return(p)
}

#' Make plot to compare predicted and reference hydrophobic thickness values
#'
#' @param df Dataframe with predicted and OPM reference values
#'
#' @return Plot object
#'
#' @examples
#' make.thk.comparison.plot(df)
#' 
make.thk.comparison.plot <- function(df) {
  for.ribbon <- seq(12, 41, by=1)
  r <- ggplot() + 
    theme_bw() + 
    background_grid() + 
    geom_abline( size = 0.35 ) + 
    geom_ribbon( aes(ymin=for.ribbon-5, ymax=for.ribbon+5, x=for.ribbon), alpha = 0.2) +
    geom_point( data = df, aes( x = opm.thk, y = best.thk, color = alpha_or_beta ), size = 0.75  ) + 
    coord_cartesian(ylim = c(15, 41), xlim =c(15,41) ) +
    scale_x_continuous( "OPM (Å)" ) +
    scale_y_continuous( "Predicted (Å)" ) +
    scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
    theme( legend.position = "none",
           text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           panel.border = element_rect( color = "black"), 
           axis.ticks = element_line( size = 0.35 ) )
  return(r)
}

orientation.frame <- function( dir, df, pdb.id, flag ) {
  
  head.df <- df
  x <- unique(df$zcoord)
  
  y <- seq( 0,90,1 )
  
  z <- seq( 0,360,5 )
  
  min.df <- data.frame()
  previous.df <- data.frame()
  energy1d.df <- data.frame()
  previous.df <- head.df[ which(head.df$azimuthal==0), ]
  #min.df <- read.table( file =  paste( dir, paste(pdb.id,"_minmatlab.dat", sep=""), sep="/"), header=T )
  min.df <- read.table( file =  paste( dir, paste(pdb.id,"_min.dat", sep=""), sep="/"), header=T )
  a <- min.df$angle[ which(min.df$angle >180 ) ]
  b <- (a - 360)
  min.df$angle[ which(min.df$angle >180 ) ] <- b
    # 
  a <- min.df$azimuthal[ which(min.df$azimuthal >180 ) ]
  b <- (a - 360)
  min.df$azimuthal[ which(min.df$azimuthal >180 ) ] <- b
  min.df$pdb.id <- pdb.id
  
 
  a <- head.df$angle[ which(head.df$angle >180 ) ]
  b <- (a - 360)
  head.df$angle[ which(head.df$angle >180 ) ] <- b
  a <- head.df$azimuthal[ which(head.df$azimuthal >180 ) ]
  b <- (a - 360)
  head.df$azimuthal[ which(head.df$azimuthal >180 ) ] <- b
  head.df$pdb.id <- pdb.id
 
  
  previous.df$pdb.id <- pdb.id
  
  a <- previous.df$angle[ which(previous.df$angle >180 ) ]
  b <- (a - 360)
  previous.df$angle[ which(previous.df$angle >180 ) ] <- b
  
  a <- previous.df$azimuthal[ which(previous.df$azimuthal >180 ) ]
  b <- (a - 360)
  previous.df$azimuthal[ which(previous.df$azimuthal >180 ) ] <- b
  
  ordered.df <- min.df[ order(min.df$total_score), ]
  write.table( ordered.df[1:100, ], file = paste( dir, paste(pdb.id,"_minordered.dat", sep=""), sep = "/") )
  minimum_depth <- ordered.df[1,1]
  minimum_azim <- ordered.df[1,3]
  energy1d.df <- head.df[ which(head.df$zcoord==minimum_depth & head.df$azimuthal==minimum_azim ), ]
  write.table( energy1d.df, file = paste( dir, paste(pdb.id,"_energy1d.dat", sep=""), sep = "/") )
 
  p1 <- make1d.energy.plots(energy1d.df)+ggtitle(paste(paste("at depth=",minimum_depth,sep=" "), paste("\n azimuthal angle=",minimum_azim,sep=" "), sep=" "))
  plot_grid( p1, ncol = 1, nrow = 1)
  a1.grid <-plot_grid( p1, ncol = 1, nrow = 1)
  #a1.grid <-plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a"), label_size = 9)
  save_plot( paste( dir, paste(pdb.id,sep="_","energy1dplot.pdf"), sep = "/"), a1.grid, units = "in", base_width = 2, base_height = 2 )
 
  ordered.df <- previous.df[ order(previous.df$total_score), ]
  write.table( ordered.df[1:100, ], file = paste( dir, paste(pdb.id,"_prevordered.dat", sep=""), sep = "/") )
  
  
  if( flag == 1 ){
    return(previous.df)
  }else{
    return(min.df)
  }
  #write.table( min.df, file = paste( dir, paste(pdb.id,"_min.dat", sep=""), sep = "/") )
 
  
}

#' Make cummulative distribution plot of the depth residuals wrt fitted line
#'
#' @param df Dataframe with predicted and OPM reference values
#'
#' @return Plot object
#'
#' @examples
#' make.cdf.depth.plot(df)
#' 
make.cdf.depth.plot <- function(df, flag) {
  
  if(flag == 0){
      s <- ggplot() + 
      theme_bw() + 
      background_grid() + 
      stat_ecdf( data = df, aes( x = abs(depth.res)), geom = "density", position = "identity", size = 0.35 ) + 
      stat_ecdf( data = df, aes( x = abs(depth.res), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) 
    }else{
      s <- ggplot() + 
        theme_bw() + 
        background_grid() + 
        stat_ecdf( data = df, aes( x = abs(depth.res_prev)), geom = "density", position = "identity", size = 0.35 ) + 
        stat_ecdf( data = df, aes( x = abs(depth.res_prev), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) 
        
   }
 s1 <- s + scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
    scale_x_continuous( "Depth Residual (Å)", limits = c(0, NA), expand = c(0,0), breaks = c(0, 5, 10, 15, 20) ) + 
    scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) + 
    theme( legend.position = "none",
           text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           panel.border = element_rect( color = "black"), 
           axis.ticks = element_line( size = 0.35 ) )
  return(s1)
}

#' Make cummulative distribution plot of the angle residuals wrt fitted line
#'
#' @param df Dataframe with predicted and OPM reference values
#'
#' @return Plot object
#'
#' @examples
#'  make.cdf.angle.plot(df)
#' 
make.cdf.angle.plot <- function(df, flag ) {
  if(flag ==0 ){
  t <- ggplot() + 
    theme_bw() + 
    background_grid() + 
    stat_ecdf( data = df, aes( x = abs(angle.res)), geom = "density", position = "identity", size = 0.35 ) + 
    stat_ecdf( data = df, aes( x = abs(angle.res), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) 
  }else{
    t <- ggplot() + 
      theme_bw() + 
      background_grid() + 
      stat_ecdf( data = df, aes( x = abs(angle.res_prev)), geom = "density", position = "identity", size = 0.35 ) + 
      stat_ecdf( data = df, aes( x = abs(angle.res_prev), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) 
    }
  
  t1<- t+ scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
    scale_x_continuous( "Angle Residual (degree ˚)", limits = c(0, NA), breaks = c(0, 5, 10, 15, 20, 25, 30, 35, 40, 45), expand = c(0,0) ) + 
    scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) + 
    theme( legend.position = "none",
           text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           panel.border = element_rect( color = "black"), 
           axis.ticks = element_line( size = 0.35 ) )
  return(t1)
}

#' Make cummulative distribution plot of the thickness residuals
#'
#' @param df Dataframe with predicted and OPM reference values
#'
#' @return Plot object
#'
#' @examples
#'  make.cdf.thk.plot(df)
#' 
make.cdf.thk.plot <- function(df) {

  df$thickness.res = (df$best.thk - df$opm.thk)
  
  u <- ggplot() + 
    theme_bw() + 
    background_grid() + 
    stat_ecdf( data = df, aes( x = abs(thickness.res)), geom = "density", position = "identity", size = 0.35 ) + 
    stat_ecdf( data = df, aes( x = abs(thickness.res), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) + 
    scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
    scale_x_continuous( "Thickness Residual (Å)", limits = c(0, NA), expand = c(0,0) ) + 
    scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) + 
    theme( legend.position = "none",
           text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           panel.border = element_rect( color = "black"), 
           axis.ticks = element_line( size = 0.35 ) )
  
  write.table( df, paste( workdir, "protein_thickness.dat", sep = "/"))
  return(u)
}

make.energy.plots <- function(df){
 
 df$zcoord = as.integer( df$zcoord )
 max_lim <- max(df$total_score[ which(abs(df$zcoord) < 35 ) ])
 
 print(max_lim)
 min_lim <- min(df$total_score[ which(abs(df$zcoord) < 35 ) ])
 
 print(min_lim)
 
 br <- c(min_lim, 0.5*(min_lim+max_lim), max_lim) 
 bestorientation.df <- compute.best.orientation( df )
  p <- ggplot() + 
    theme_bw() + 
    geom_raster( data = df, aes( x = angle, y = zcoord, fill = total_score) ) +
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(0,45,90,135,180) )+
    scale_y_continuous( "Depth (A)", expand = c(0,0), limits = c(-30, 30) ) + 
    scale_fill_viridis( "Energy (REU)", limits=c( min_lim, max_lim ), breaks = c( min_lim, max_lim ) ) + 
    geom_point( data = bestorientation.df, aes( x = best.tilt, y = best.depth ), color = "white" ) + 
    geom_point( data = bestorientation.df, aes( x = best.tilt, y = best.depth ), color = "white", shape = 1, size = 1.0 ) +
    geom_point( data = opm.orientation[ which( opm.orientation$pdb == unique( df$pdb.id ) ), ], aes( x = angle, y = depth ), color = "red")+
    geom_point( data = opm.orientation[ which( opm.orientation$pdb == unique( df$pdb.id ) ), ], aes( x = angle, y = depth ), color = "red", shape = 2, size = 2.0 ) + 
    theme( legend.position = "bottom",
           text = element_text( size = 7, color = "black"), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35, color = "black"), 
           panel.border = element_rect(),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-2,-7))
 
  return(p)
}

make1d.energy.plots <- function(df){
  
  df$zcoord = as.integer( df$zcoord )
  max_lim = max(df$total_score) + 10
  min_lim = min(df$total_score) - 10
  
  p <- ggplot( data = df, aes(x=angle, y=total_score) ) + 
    geom_line( color="grey") +
    geom_point(shape=21, color="black", fill="#69b3a2", size=2) +
    theme_bw() +
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(0,45,90,135,180) )+
    scale_y_continuous( "Energy (REU)", expand = c(0,0), limits = c(min_lim, max_lim)) + 
    theme( legend.position = "bottom",
           text = element_text( size = 7, color = "black"), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35, color = "black"), 
           panel.border = element_rect(),
           legend.key.size = unit(0.5, "lines"), 
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,-7,-2,-7))
  ## scale_fill_viridis( "Energy (REU)", limits=c( minlim, maxlim ), breaks = c( minlim, maxlim ) )
  return(p)
}



workdir <- "/path-to/protein-tilt-angle"
workdir_hl <- "/path-to/hydrophobic-length"

experimental.data <- "experiment_test3.txt" 
#experimental file in the workdir

opm.orientation <- read.table( paste( workdir, experimental.data, sep = "/"), header = T )
main.df <- read.mappings( workdir, experimental.data )
main.df2 <- read.mappings.hl( workdir_hl, experimental.data )

residual.df <- compute.and.output.residuals( main.df )

p1<-make.tilt.angle.comparison.plot(main.df, 0) + ggtitle("Minimize over azimuthal angle")
p2<-make.tilt.angle.comparison.plot(main.df, 1) + ggtitle("At azimuthal angle = 0" )

a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
save_plot( paste( workdir, paste("tilt",sep="_","compare.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )

p1<-make.depth.comparison.plot(main.df, 0) + ggtitle("Minimize over azimuthal angle")
p2<-make.depth.comparison.plot(main.df, 1) + ggtitle("At azimuthal angle = 0" )

a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
save_plot( paste( workdir, paste("depth",sep="_","compare.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )

p1<-make.cdf.angle.plot(residual.df, 0) + ggtitle("Minimize over azimuthal angle")
p2<-make.cdf.angle.plot(residual.df, 1) + ggtitle("At azimuthal angle = 0" )

a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
save_plot( paste( workdir, paste("tilt",sep="_","freqcompare.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )

p1<-make.cdf.depth.plot(residual.df, 0) + ggtitle("Minimize over azimuthal angle")
p2<-make.cdf.depth.plot(residual.df, 1) + ggtitle("At azimuthal angle = 0" )

a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
save_plot( paste( workdir, paste("depth",sep="_","freqcompare.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )

p1<-make.cdf.thk.plot(main.df2) + ggtitle("CDF plot" )
p2<-make.thk.comparison.plot(main.df2) + ggtitle("Scatter plot" )

a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
save_plot( paste( workdir, paste("length",sep="_","compare_wrtexp_adjusted.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )

