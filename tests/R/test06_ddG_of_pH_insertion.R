#' Analyze dG of peptide insertion results for different pH values (Test #6). 
# @author: Rituparna Samanta (rituparna@utexas.edu)  and Rebecca Alford <ralford3[at]jhu.edu> 
# @brief: Computes analysis for peptide insertion calculations, including:
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
#' @return out.df Dataframe with best tilt angle and depth position 
#'
#' @examples
#' compute.best.membrane.orientation(df)
#'
compute.best.membrane.orientation <- function(df) {
  
  # Make a dataframe with only the membrane positions
  memb.df <- df[ which( abs(df$zcoord) <= 15 ), ]
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
  print(best.soluble.tilt)
  print(best.soluble.depth)
  # Create a dataframe that stores all of this information and return it
  out.df <- data.frame( best.soluble.tilt = best.soluble.tilt,
                        best.soluble.depth = best.soluble.depth,
                        best.soluble.score = best.soluble.score )
  #best.soluble.azim = best.soluble.azim,
  return(out.df)
}

#' Read the energy landscapes for different peptides
#'
#' @param path to ddG-of-pH-insertion directory
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
    data.file <- list.files(path = data.dirs[d], pattern = "*landscape.sh")
    # Get the file metadata
    file.info <- strsplit(data.file, "_")[1]
    pdb.id <- unlist(file.info)[1]
    print(pdb.id)
    dir.info <- strsplit(data.dirs[d], "_")
    pH.value <- unlist(dir.info)[length(unlist(dir.info))]
    print(pH.value)
    
    file.name <- paste(data.dirs[d],paste(pdb.id,"combined.dat", sep="_"), sep = "/")
    sub.df <- read.table( file.name, header = T )
    
    
    dG.pH4 <- 0
    dG.pH8 <- 0
    
    sub_prev.df <- orientation.frame( sub.df, 1, pdb.id, data.dirs[d], pH.value )
    sub_min.df <- orientation.frame( sub.df, 0, pdb.id, data.dirs[d], pH.value )
  
  
     if ( pH.value == 4 ) {
       membrane.orientation <- compute.best.membrane.orientation(sub_min.df)
       dG.pH4 <- membrane.orientation$best.memb.score
       pH4.df <- data.frame( pdb.id = pdb.id, pH.value = pH.value, dG = dG.pH4, tilt = membrane.orientation$best.memb.tilt, depth = membrane.orientation$best.memb.depth )
       #making a list of pdb.id,pH.value, score at best orientation in membrane
       main.df <- rbind( main.df, pH4.df )
    } else {
      soluble.orientation <- compute.best.soluble.orientation(sub_min.df)
      dG.pH8 <- soluble.orientation$best.soluble.score
      pH8.df <- data.frame( pdb.id = pdb.id, pH.value = pH.value, dG = dG.pH8, tilt = soluble.orientation$best.soluble.tilt, depth = soluble.orientation$best.soluble.depth )
        #making a list of pdb.id,pH.value, score at best orientation outside membrane
      main.df <- rbind( main.df, pH8.df )
    }
    
     
    p1 <- make.energy.plots(sub_min.df)
    p2 <- make.energy.plots(sub_prev.df)
    plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
    a1.grid <-plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
    save_plot( paste( data.dirs[d], paste(pdb.id,sep="_","full.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )
    
    
     if ( pH.value == 4 ) {
       membrane.orientation <- compute.best.membrane.orientation( sub_prev.df )
       dG.pH4 <- membrane.orientation$best.memb.score
       pH4.df <- data.frame( pdb.id = pdb.id, pH.value = pH.value, dG = dG.pH4, tilt = membrane.orientation$best.memb.tilt, depth = membrane.orientation$best.memb.depth )
 
       main_prev.df <- rbind( main_prev.df, pH4.df )
     } else {
       soluble.orientation <- compute.best.soluble.orientation(sub_prev.df)
       dG.pH8 <- soluble.orientation$best.soluble.score
       pH8.df <- data.frame( pdb.id = pdb.id, pH.value = pH.value, dG = dG.pH8, tilt = soluble.orientation$best.soluble.tilt, depth = soluble.orientation$best.soluble.depth )
    
       main_prev.df <- rbind( main_prev.df, pH8.df )
     }
     
   
 
   write.table( main_prev.df, file = paste( workdir, "previous_dG_pH.txt", sep = "/") )
   write.table( main.df, file = paste( workdir, "min_dG_pH.txt", sep = "/") )
  
  } 
  
    data_prev.df <- data.frame()
    data.df <- data.frame()
    
    main_prev.df <- data.frame(read.table(file = paste( workdir, "previous_dG_pH.txt", sep = "/"), header=T))
    main.df <- data.frame(read.table(file = paste( workdir, "min_dG_pH.txt", sep = "/"), header=T))
    seq.tags <- unique(main_prev.df$pdb.id)
   #### all pdb
    for (seq in seq.tags) {
      
       print(seq)
       dG.pH4 <- main.df$dG[ which( main.df$pH.value == 4 & main.df$pdb.id == seq )]
       dG.pH8 <- main.df$dG[ which( main.df$pH.value == 8 & main.df$pdb.id == seq )]
       ddG.of.pH.shift = round( dG.pH4 - dG.pH8, 3)
       variant.num <- unlist(strsplit(as.character(seq), "-|\\s")[1])[2]
       exp.ddG <- experimental.data$ddG[ which( experimental.data$Sequence == seq ) ]
    ###   #exp.ddG = dG_pH4 - dG_pH8
       sub.df <- data.frame( seq = seq, var.no = variant.num, pred.ddG = ddG.of.pH.shift, exp.ddG = exp.ddG )
       data.df <- rbind( data.df, sub.df )
       
       dG1.pH4 <- main_prev.df$dG[ which( main_prev.df$pH.value == 4 & main_prev.df$pdb.id == seq )]
       dG1.pH8 <- main_prev.df$dG[ which( main_prev.df$pH.value == 8 & main_prev.df$pdb.id == seq )]
       ddG1.of.pH.shift = round( dG1.pH4 - dG1.pH8, 3)
       variant.num <- unlist(strsplit(as.character(seq), "-|\\s")[1])[2]
     ## #     print(variant.num)
       exp.ddG <- experimental.data$ddG[ which( experimental.data$Sequence == seq ) ]
     ##  #exp.ddG = dG_pH4 - dG_pH8
       sub.df <- data.frame( seq = seq, var.no = variant.num, pred.ddG = ddG1.of.pH.shift, exp.ddG = exp.ddG )
       data_prev.df <- rbind( data_prev.df, sub.df )
       
    }
    
    linearmod <-lm(pred.ddG~exp.ddG, data.df)
    print("details for minimum")
    print(linearmod)
    summary(linearmod)
   
    linearmod_prev <-lm(pred.ddG~exp.ddG, data_prev.df)
    print("details for minimum")
    print(linearmod_prev)
    summary(linearmod_prev)
    
    write.table( data.df, file = paste( workdir, "processed_ddG_pH.txt", sep = "/") )
    write.table( data_prev.df, file = paste( workdir, "processed_ddG_pH_prev.txt", sep = "/") )
 
    p1 <- make.pH.dG.ins.correlation.plot( data.df, as.numeric(linearmod$coefficients[2]), as.numeric(linearmod$coefficients[1]) )
  ##  #r2-0.085; r2 after alignment = 0.03309 
    p2 <- make.pH.dG.ins.correlation.plot( data_prev.df, as.numeric(linearmod_prev$coefficients[2]), as.numeric(linearmod_prev$coefficients[1]) )
    ###r2 after alignment = 0.03826
    plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
    a1.grid <-plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
    save_plot( paste( workdir, paste("processed",sep="_","full.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 2 )
    
}

#' Extract the dataframe of energy landscape minimized over all azimuthal angle and at azimuthal angle =0 
#'
#' @param directory containing the *combined.dat, *min.dat; datfarame with energy at depth and angle; 
#' flag about the options, pdb
#'
#' @return filtered data frame
#'
#' @examples
#' orientation.frame(df, 0, 'pHLIP-v1','dir', 4)
#' 
orientation.frame <- function( df, flag, pdb.id, dir, pH.value ) {
  
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
make.pH.dG.ins.correlation.plot <- function( pH.df, pH.slope, pH.intercept ) {
 # temp <- locator(1) # On the chart, click where you would like the text to appear
  max_lim = max(pH.df$pred.ddG)
  min_lim = min(pH.df$pred.ddG)
  p <- ggplot( data = pH.df, aes( x = exp.ddG, y = pred.ddG, label = var.no ) ) +
    theme_bw() + 
    background_grid() +
    geom_abline( slope = pH.slope, intercept = pH.intercept, color = "gray50", size = 0.35 ) + 
    geom_point( size = 0.75 ) +
    geom_text_repel( size = 2.5, force = 3 ) +
  #  legend("bottomright",legend= paste("Slope = ",pH.slope), bty ="n", pch=NA)+
    scale_x_continuous( "Experiment (kcal/mol)", limits = c(0, 3.0), expand = c(0,0) ) + 
    scale_y_continuous( "Predicted (REU)" , limits = c(min_lim, max_lim), expand = c(0,0) ) + 
    theme( text = element_text( size = 7, color = "black" ), 
           axis.text = element_text( size = 7, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( size = 0.35, color = "black" ), 
           panel.border = element_rect( color = "black" ))
  return(p)
 
}

#' Make energy landscape plots
#'
#' @param df Dataframe with energy landscape
#'
#' @return Plot Object
#'
#' @examples
#' make.energy.plots(df)
#' 
make.energy.plots <- function(df){
    df$zcoord = as.integer( df$zcoord )
    p <- ggplot() + 
    theme_bw() + 
    geom_raster( data = df, aes( x = angle, y = zcoord, fill = total_score ) ) +
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) + 
    scale_y_continuous( "Depth (Å)", expand = c(0,0), limits = c(-45, 45) ) + 
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



workdir <- "/path-to-/ddG-of-pH-insertion"
experimental.data <- read.table( paste( workdir, "experiment_test6.dat", sep = "/"), header = T )

read.mappings(workdir)

