#' Analyze dG of peptide insertion results (Test #5)
#' @author: Rebecca Alford <ralford3[at]jhu.edu> and Rituparna Samanta <rsamant2[at]jhu.edu>
#' @brief: Computes analysis for peptide insertion calculations, including:
#'    - Correlation between predicted and experimentally measured value

library(nnet)
library(viridis)
library(ggplot2)
library(ggrepel)
library(cowplot)

#' Compute the lowest energy orientation inside the membrane
#'
#' @param df Dataframe with mapping between peptide orientations and energies
#'
#' @return Dataframe with best tilt angle and depth position 
#'
#' @examples
#' compute.best.membrane.orientation(df)
#'
compute.best.membrane.orientation <- function(df) {
  
  # Make a dataframe with only the membrane positions
  memb.df <- df[ which( abs(df$zcoord) <= 12 ), ]
  ordered.memb.df <- memb.df[ order(memb.df$total_score), ]
  head <- head(ordered.memb.df,1)
  best.memb.tilt <- head$angle
  #best.memb.azim <- head$azimuthal
  best.memb.depth <- head$zcoord
  best.memb.score <- head$total_score
  print('membrane')
  print(best.memb.tilt)
  print(best.memb.depth)
  print(best.memb.score)
  # Create a dataframe that stores all of this information and return it
  out.df <- data.frame( best.memb.tilt = best.memb.tilt, 
                        best.memb.depth = best.memb.depth,
                        best.memb.score = best.memb.score )
  # best.memb.azim = best.memb.azim, 
  return(out.df)
}

#' Compute the lowest energy orientation outside the membrane
#'
#' @param df Dataframe with mapping between peptide orientations and energies
#'
#' @return Dataframe with best tilt angle and depth position 
#'
#' @examples
#' compute.best.soluble.orientation(df)
#'
compute.best.soluble.orientation <- function(df) {
  
  # Make a dataframe with only the membrane positions
  memb.df <- df[ which( abs(df$zcoord) > 30 ), ]
  ordered.memb.df <- memb.df[ order(memb.df$total_score), ]
  head <- head(ordered.memb.df,1)
  best.soluble.tilt <- head$angle
  #best.soluble.azim <- head$azimuthal
  best.soluble.depth <- head$zcoord
  best.soluble.score <- head$total_score
  print('soluble')
  print(best.soluble.score)
  # Create a dataframe that stores all of this information and return it
  out.df <- data.frame( best.soluble.tilt = best.soluble.tilt,
                        best.soluble.depth = best.soluble.depth,
                        best.soluble.score = best.soluble.score )
  #best.soluble.azim = best.soluble.azim,
  return(out.df)
}

#' Read the energy landscapes for different peptides
#'
#' @param path to ddG-of-insertion directory
#'
#' @return  processed_data.txt files with ddG of insertion for different conditions. 
#'
#' @examples
#' read.mappings('dir')
#' 

read.mappings <- function( dir ) {
  
  main.df <- data.frame()
  main_prev.df <- data.frame()
  
  #input : dir is the work directory
  data.dirs <- list.dirs(path = dir)
  for(d in 2:(length(data.dirs))) {
      #
      # Get the file
      data.file <- list.files(path = data.dirs[d], pattern = "*combined.dat")
      # Get the file metadata
      file.info <- strsplit(data.file, "_")[1]
      pdb.id <- unlist(file.info)[1]
      print(pdb.id)
      dir.info <- strsplit(data.dirs[d], "_")
     
      sub.df <- data.frame()
       
      for( index in 1:length(data.file) ){
          file.name <- paste(data.dirs[d], data.file[index], sep = "/")
          df <- read.table( file.name, header = T )
          sub.df <- rbind(sub.df,df)
      }
      dG.in <- 0
      dG.out <- 0
      sub_prev.df <- orientation.frame( data.dirs[d], sub.df, 1, pdb.id )
      sub_min.df <- orientation.frame( data.dirs[d], sub.df, 0, pdb.id )
      
     
       membrane.orientation <- compute.best.membrane.orientation(sub_min.df)
       dG.in <- membrane.orientation$best.memb.score
       in.df <- data.frame( pdb.id = pdb.id, location = "membrane", dG = dG.in, tilt = membrane.orientation$best.memb.tilt, depth = membrane.orientation$best.memb.depth )
       #making a list of pdb.id,pH.value, score at best orientation in membrane
       main.df <- rbind( main.df, in.df )
    
       soluble.orientation <- compute.best.soluble.orientation(sub_min.df)
       dG.out <- soluble.orientation$best.soluble.score
       out.df <- data.frame( pdb.id = pdb.id, location = "soluble", dG = dG.out, tilt = soluble.orientation$best.soluble.tilt, depth = soluble.orientation$best.soluble.depth )
       #making a list of pdb.id,pH.value, score at best orientation outside membrane
       main.df <- rbind( main.df, out.df )
    
     
      p1 <- make.energy.plots(sub_min.df)
      p2 <- make.energy.plots(sub_prev.df)
      plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
      a1.grid <-plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
      save_plot( paste( data.dirs[d], paste(pdb.id,sep="_","full.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )
      
    
    
      membrane.orientation <- compute.best.membrane.orientation( sub_prev.df )
      dG.in <- membrane.orientation$best.memb.score
      in.df <- data.frame( pdb.id = pdb.id, location = "membrane", dG = dG.in, tilt = membrane.orientation$best.memb.tilt, depth = membrane.orientation$best.memb.depth )
      #making a list of pdb.id,pH.value, score at best orientation in membrane
      main_prev.df <- rbind( main_prev.df, in.df )
    
      soluble.orientation <- compute.best.soluble.orientation(sub_prev.df)
      dG.out <- soluble.orientation$best.soluble.score
      out.df <- data.frame( pdb.id = pdb.id, location = "soluble", dG = dG.out, tilt = soluble.orientation$best.soluble.tilt, depth = soluble.orientation$best.soluble.depth )
      #making a list of pdb.id,pH.value, score at best orientation outside membrane
      main_prev.df <- rbind( main_prev.df, out.df )
  } 
    
  
#  print('writing')
  
  write.table( main_prev.df, file = paste( workdir, "previous_dG.txt", sep = "/") )
  write.table( main.df, file = paste( workdir, "min_dG.txt", sep = "/") )
  
  # Calculate the ddG of insertion values
   data_prev.df <- data.frame()
   data.df <- data.frame()
   seq.tags <- unique(main.df$pdb.id)
  # all pdb
   for (seq in seq.tags) {
     
     print(seq)
      dG.in <- main.df$dG[ which( main.df$location == "membrane" & main.df$pdb.id == seq )]
      dG.out <- main.df$dG[ which( main.df$location == "soluble" & main.df$pdb.id == seq )]
      ddG.of.shift = round( dG.in - dG.out, 3)
     
      exp.ddG <- experimental.data$ddG_insert[ which( experimental.data$Sequence == seq ) ]
      #exp.ddG = dG_pH4 - dG_pH8
      sub.df <- data.frame( seq = seq, pred.ddG = ddG.of.shift, exp.ddG = exp.ddG )
      data.df <- rbind( data.df, sub.df )
      
      dG1.in <- main_prev.df$dG[ which( main_prev.df$location == "membrane" & main_prev.df$pdb.id == seq )]
      dG1.out <- main_prev.df$dG[ which( main_prev.df$location == "soluble" & main_prev.df$pdb.id == seq )]
      ddG1.of.shift = round( dG1.in - dG1.out, 3)
      #     print(variant.num)
      exp.ddG <- experimental.data$ddG_insert[ which( experimental.data$Sequence == seq ) ]
      #exp.ddG = dG_pH4 - dG_pH8
      sub.df <- data.frame( seq = seq, pred.ddG = ddG1.of.shift, exp.ddG = exp.ddG )
      data_prev.df <- rbind( data_prev.df, sub.df )
      
   }
   write.table( data.df, file = paste( workdir, "processed_ddG.txt", sep = "/") )
   write.table( data_prev.df, file = paste( workdir, "processed_ddG_prev.txt", sep = "/") )
 
  
}
#' Extract the dataframe of energy landscape minimized over all azimuthal angle and at azimuthal angle =0 
#'
#' @param directory containing the *combined.dat, *min.dat; datfarame with energy at depth and angle; 
#' flag about the options, pdb
#'
#' @return filtered data frame
#'
#' @examples
#' orientation.frame('dir', df, 0, 'GL5')
#' 
orientation.frame <- function( dir, df, flag, pdb.id ) {
  
  x <- unique(df$zcoord)
  # seq( min(df$zcoord), max(df$zcoord), 1 )
  y <- seq( 0,360,1 )
  
  z <- seq( 0,360,5 )
  
  
  if( flag == 0 ){
    min.df <- data.frame()
    
    min.df <- read.table( file = paste(dir, paste(pdb.id,"_min.dat", sep="") ,sep="/"), header=T )
    
    a <- min.df$angle[ which(min.df$angle > 180 ) ]
    b <- (a - 360)
    min.df$angle[ which(min.df$angle > 180 ) ] <- b
    
    a <- min.df$azimuthal[ which(min.df$azimuthal > 180 ) ]
    b <- (a - 360)
    min.df$azimuthal[ which(min.df$azimuthal > 180 ) ] <- b
    
    return( min.df )
  } else if( flag ==1 ){
    
    previous.df <- data.frame()
    previous.df <- df[ which(df$azimuthal==0), ]
    
    a <- previous.df$angle[ which(previous.df$angle > 180 ) ]
    b <- (a - 360)
    previous.df$angle[ which(previous.df$angle > 180 ) ] <- b
    
    a <- previous.df$azimuthal[ which(previous.df$azimuthal > 180 ) ]
    b <- (a - 360)
    previous.df$azimuthal[ which(previous.df$azimuthal > 180 ) ] <- b
    
    return( previous.df )
  } else{
    print( 'wrong option' )
  }
  
}

#' Make correlation plot between predicted and experimentally measured values
#'
#' @param df Dataframe with predicted and experimentally measured values
#'
#' @return Plot Object
#'
#' @examples
#' make.pH.dG.ins.correlation.plot(df)
#' 
make.dG.ins.correlation.plot <- function( df, slope, intercept ) {
 # temp <- locator(1) # On the chart, click where you would like the text to appear
  linearmod <-lm(pred.ddG~exp.ddG, df)
  print(linearmod)
 
 
  p <- ggplot( data = df, aes( x = exp.ddG, y = pred.ddG, label = seq ) ) +
    theme_bw() + 
    background_grid() +
    geom_vline( xintercept = 0, color = "gray50", linetype = "dashed", size = 0.35 ) + 
    geom_abline( slope = slope, intercept = intercept, color = "gray50", size = 0.35 ) + 
    geom_point( size = 0.75 ) +
   # geom_errorbarh( aes( xmax = exp+error, xmin = exp-error ) ) +
    geom_text_repel( size = 2.5, force = 3 ) +
  #  legend("bottomright",legend= paste("Slope = ",pH.slope), bty ="n", pch=NA)+
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(-2.5, 2.5), expand = c(0,0) ) + 
    scale_y_continuous( "Predicted (REU)" , limits = c(-20, -10), expand = c(0,0) ) + 
    theme( text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35, color = "black" ), 
           panel.border = element_rect( color = "black" ))
  return(p)
  #how to put a text in the figure? 
  # geom_abline( slope = 1.338, intercept = -11.15, color = "gray50", size = 0.35 ) + 
  # geom_text_repel( size = 2.5, force = 3 ) + 
  #, limits = c(-20, -5), expand = c(0,0)
}

make.energy.plots <- function(df){
  p <- ggplot() + 
    theme_bw() + 
    geom_raster( data = df, aes( x = angle, y = zcoord, fill = total_score ) ) +
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) + 
    scale_y_continuous( "Depth (A)", expand = c(0,0), limits = c(-45, 45) ) + 
    scale_fill_viridis( "Energy\n(REU)" ) + 
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



workdir <- "/path-to/ddG-of-insertion"
experimental.data <- read.table( paste( workdir, "experiment_test5.dat", sep = "/"), header = T )

read.mappings(workdir)
data.df <- read.table( paste( workdir,"processed_ddG.txt", sep = "/"), header = T )
data_prev.df <- read.table( paste( workdir, "processed_ddG_prev.txt", sep = "/"), header = T )

#p1 <- make.dG.ins.correlation.plot( data.df, 2.41, -16.01 ) + ggtitle(label = "Minimum over azimuthal angle")
##R2 = 0.9911
p1 <- make.dG.ins.correlation.plot( data.df, 2.499, -16.876 ) + ggtitle(label = "Minimum over azimuthal angle")
#R2 = 0.9911
#p2 <- make.dG.ins.correlation.plot( data_prev.df, 2.407, -15.783 ) + ggtitle(label = "At azimuthal angle = 0")
p2 <- make.dG.ins.correlation.plot( data_prev.df, 2.472, -16.703 ) + ggtitle(label = "At azimuthal angle = 0")

a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)

save_plot( paste( workdir, paste("comparing_exp_simulation_test5",sep="_","full.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )
