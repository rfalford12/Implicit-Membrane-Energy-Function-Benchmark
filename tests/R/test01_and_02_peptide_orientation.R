#' @author RituparnaSamanta <rsamant2[at]jhu.edu> and Rebecca Alford <ralford3[at]jhu.edu> 
# @brief: Compare experimental and predicted tilt angle (minimum over rotation angle) and depth (zcoord)

library(nnet)
library(viridis)
library(ggplot2)
library(cowplot)

#' Read the *_combined.dat file
#'
#' @param Workdir name and the folder name in which the _combined.dat are kept
#'
#' @return Dataframe with tilt angle, depth position, total_score
#'
#' @examples
#' read.mapping('dir', 'combined_folder')

read.mapping <- function( dir, test.name ) {
  
  # Create a dataframe with the energy landscapes
  main.df <- data.frame()
  
  # List out all data files in the test category
  test.path = paste( dir, test.name, sep = "/" )
  
  data.files <- list.files(path = test.path, pattern = "*combined.dat")
 
  for(f in 1:length(data.files)) {
    
    # Grab the pdb identifier
    pdb.id <- strsplit(data.files[f], "_combined")[1]
    pdb.id <- unlist(pdb.id)[1]
    print(pdb.id)
    # Read data from landscape file
    file.name <- paste(test.path, data.files[f], sep = "/")
    df<- data.frame(read.table( file.name, header = T ))
    
    print(data.files[f])
    #df <-data.frame(read.table(file=file.name,header=T))
    
    
   
    a <- df$angle[ which(df$angle >180 ) ]
    b <- (a - 360)
    df$angle[ which(df$angle >180 ) ] <- b
    
    a <- df$azimuthal[ which(df$azimuthal >180 ) ]
    b <- (a - 360)
    df$azimuthal[ which(df$azimuthal >180 ) ] <- b
    
    df$pdb.id <- pdb.id
    
    main.df <- rbind( main.df, df )
  }
  return(main.df)
}

#' extracts the subset of the energy landscape of a particular pdb.id
#'
#' @param df, pdb.id and the folder name in which the _combined.dat are kept
#'
#' @return Dataframe extracts energy landscape at depth=+/-15A and at azimuthal angle = 0
#'
#' @examples
#' orientation.frame(df, '1a11', 'combined_folder')

orientation.frame <- function( df, pdb.id, dir ) {
  
  head.df <- data.frame()
  head.df <- df[ which( df$pdb.id == pdb.id ), ]
  
  #x <- seq( head.df$zcoord[1], head.df$zcoord[nrow(head.df)], 1 )
  x <- unique( head.df$zcoord )
  
  
  y <- seq( 0,360,1 )
  
  
  z <- seq( 0,360,5 )
  
  
  zat15.df <- data.frame()
  zatm15.df <- data.frame()
  zat15.df <- head.df[ which( head.df$zcoord == 15 ), ]
  zatm15.df <- head.df[ which( head.df$zcoord == -15 ), ]
  
  min.df <- data.frame()
  previous.df <- data.frame()
  previous.df <- head.df[ which(head.df$azimuthal==0), ]
  
  
  min.df <- read.table( file = paste(dir, paste(pdb.id,"_min.dat", sep="") ,sep="/"), header=T )

  #to avoid two repeating results for normal=0
  a <- min.df$angle[ which(min.df$angle >180 ) ]
  b <- (a - 360)
  min.df$angle[ which(min.df$angle >180 ) ] <- b

  a <- min.df$azimuthal[ which(min.df$azimuthal >180 ) ]
  b <- (a - 360)
  min.df$azimuthal[ which(min.df$azimuthal >180 ) ] <- b
  min.df$pdb.id <- pdb.id

   
  ordered.df <- min.df[ order(min.df$total_score), ]
  #if( flag == 0 ){
  write.table( min.df, file = paste( dir, paste(pdb.id,"_min.dat", sep=""), sep = "/") )
  write.table( min.df[ ,1:4 ], file = paste( dir, paste(pdb.id,"_min.csv", sep=""), sep = "/"), sep =",",na="NA",row.names=FALSE, col.names=FALSE)
  write.table( ordered.df[ which( ordered.df$total_score <= (ordered.df[1,4]+0.20) ), ], file = paste( dir, paste(pdb.id,"_minordered.dat", sep=""), sep = "/") )
  
  min_zcoord = ordered.df[1,1]
  min_angle = ordered.df[1,2]
  if(min_zcoord<0){ 
    min_limit = -40
    max_limit = 0
  }else
  {
    max_limit = 40
    min_limit = 0
  }
  if(min_angle<0){
    min_angle_lim = -90
    max_angle_lim = 0
  }else{
    min_angle_lim = 0
    max_angle_lim = 90
    
  }
  ordered.df <- ordered.df[ which( ordered.df$total_score <= (ordered.df[1,4]+0.10) & ( ordered.df$angle>=min_angle_lim & ordered.df$angle<=max_angle_lim ) & (ordered.df$zcoord>=min_limit & ordered.df$zcoord<=max_limit) ), ]
  write.table( ordered.df, file = paste( dir, paste(pdb.id,"_minorderedfirstphase.dat", sep=""), sep = "/") )
  print("min is:")
  print(min(ordered.df$angle))
  print("max is:")
  print(max(ordered.df$angle))
  
 #_previous is at azimuthal angle = 0 
  ordered.df <- previous.df[ order(previous.df$total_score), ]
  write.table( previous.df, file = paste( dir, paste(pdb.id,"_previous.dat", sep=""), sep = "/") )
  write.table( previous.df[ ,1:4 ], file = paste( dir, paste(pdb.id,"_previous.csv", sep=""), sep = "/"), sep=",",na="NA",row.names=FALSE, col.names=FALSE)
  write.table( ordered.df[ which( ordered.df$total_score <= (ordered.df[1,4]+2.0) ), ], file = paste( dir, paste(pdb.id,"_prevordered.dat", sep=""), sep = "/") )
  
  write.table( zat15.df, file = paste( dir, paste(pdb.id,"_at15.dat", sep=""), sep = "/") )
  write.table( zat15.df[ ,1:4 ], file = paste( dir, paste(pdb.id,"_at15.csv", sep=""), sep = "/"), sep=",",na="NA",row.names=FALSE, col.names=FALSE)
  
  write.table( zatm15.df, file = paste( dir, paste(pdb.id,"_atm15.dat", sep=""), sep = "/") )
  write.table( zatm15.df[ ,1:4 ], file = paste( dir, paste(pdb.id,"_atm15.csv", sep=""), sep = "/"), sep=",", na="NA",row.names=FALSE, col.names=FALSE)
  
}

#' extracts the best orientation when energy landscape 
#'
#' @param df, flag=2,3 when the depths are constant; 0,1 when all variables are considered. 
#'
#' @return based on the flag, it calculates the tilt, angle and depth at which there is minimum score.  
#'
#' @examples
#' compute.best.orientation(df, 0)

compute.best.orientation <- function( df, flag ) { 
  df.out = data.frame()
  ordered.df <- df[ order(df$total_score), ]
  head <- head(ordered.df,1)
  
  if( flag == 0 || flag == 1){
   df.out = data.frame( best.x = ordered.df$angle[ which( ordered.df$total_score == head$total_score )], best.y = ordered.df$zcoord[ which( ordered.df$total_score == head$total_score )], best.z = ordered.df$azimuthal[ which( ordered.df$total_score == head$total_score )])
  } else if( flag == 2 || flag ==3 ){
   df.out = data.frame( best.x = ordered.df$angle[ which( ordered.df$total_score == head$total_score )], best.y = ordered.df$azimuthal[ which( ordered.df$total_score == head$total_score )])
  } 
  
  #print(df.out)
  return(df.out)
}

#' plot the energy landscape 
#'
#' @param dir, flag=2,3 when the depths are constant; 0,1 when all variables are considered, pdb.id 
#'
#' @return based on the flag, it plots the energy landcsapes  
#'
#' @examples
#' plot.orientation.map(dir, 0, '1a11')
#' 
plot.orientation.map <- function( dir, flag, pdb.id ) {
  
  
  if( flag == 0 || flag == 1){
  
    df_min<- data.frame(read.table( file = paste( dir, paste(pdb.id,"_min.dat", sep=""), sep = "/"), header=T ))
    df_prev <- data.frame(read.table( file = paste( dir, paste(pdb.id,"_previous.dat", sep=""), sep = "/"), header=T ))
    
    maxlim_min <- max( df_min$total_score )
    minlim_min <- min( df_min$total_score )
    maxlim_prev <- max( df_prev$total_score )
    minlim_prev <- min( df_prev$total_score )
    
    maxlim <- max(maxlim_min, maxlim_prev)
    minlim <- min(minlim_min, minlim_prev)
    
    if(flag==0){
      df <- df_min
      
    }else{
      df <- df_prev
    }
    #  return( min.df )
  } else if( flag == 2 || flag == 3){
    df_z15 <- data.frame(read.table( file = paste( dir, paste(pdb.id,"_at15.dat", sep=""), sep = "/"), header=T ))
    df_zm15 <- data.frame( read.table( file = paste( dir, paste(pdb.id,"_atm15.dat", sep=""), sep = "/"), header=T ))
    
    maxlim_z15 <- max( df_z15$total_score )
    minlim_z15 <- min( df_z15$total_score )
    maxlim_zm15 <- max( df_zm15$total_score )
    minlim_zm15 <- min( df_zm15$total_score )
    
    maxlim <- max( maxlim_z15, maxlim_zm15 )
    minlim <- min( minlim_z15, minlim_zm15 )
    
    if( flag == 2){
      df <- df_z15
    }else{
      df <- df_zm15
    }
    
  }else{ 
    print( "wrong choice" )
  }
  
  
 

   # Get coordinates with the best orientation
  bestorientation.df <- compute.best.orientation( df, flag )
  #print( bestorientation.df )
  df$zcoord = as.integer( df$zcoord )
 if(flag==0 || flag==1)
  {p <- ggplot() + 
    theme_bw() + 
    geom_raster( data = df, aes( x = angle, y = zcoord, fill = total_score ) ) +
    geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white" ) + 
    geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white", shape = 1, size = 0.5 ) + 
    geom_point( data = experimental.data[ which( experimental.data$pdb == pdb.id ), ], aes( x = tilt, y = z_coord ), color = "red" ) + 
    geom_point( data = experimental.data[ which( experimental.data$pdb == pdb.id ), ], aes( x = tilt, y = z_coord ), color = "red", shape = 2, size = 2 ) + 
    scale_x_continuous( "Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) + 
    scale_y_continuous( "zcoord (Å)", expand = c(0,0), limits = c(-45, 45) ) + 
    scale_fill_viridis( "Energy\n(REU)" , limits=c( minlim, maxlim)) +
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
else{
     df1 <- df[order(df$angle, df$azimuthal),]
     df=df1
    p <- ggplot( ) + 
      theme_bw() + 
      geom_raster( data = df, aes( x = angle, y = azimuthal, fill = total_score ) ) +
      geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white" ) + 
      geom_point( data = bestorientation.df, aes( x = best.x, y = best.y ), color = "white", shape = 1, size = 0.5 ) + 
      scale_x_continuous( "Normal Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) + 
      scale_y_continuous( "Azimuthal Angle (degrees)", expand = c(0,0), breaks = c(-180, -90, 0, 90, 180) ) + 
      scale_fill_viridis( "Energy\n(REU)" , limits=c( minlim, maxlim)) +
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
  

}
workdir <- "/path-to-workdir/"
experimental.data <- read.table( paste( workdir, "experimental_tilt.dat", sep = "/"), header = T )

#pass the name of the 'combined_folder' which contains all pdb.id_combined.dat and pdb.id_min.dat
df.tm.natives <- read.mapping( workdir, 'git_test1and2' )



pdb_full <- c("1a11", "1mp6", "2nr1", "1pje", "WALP23", "polyA-cappedW", "polyA-cappedY", "2mag", "LK-peptide-n6", "1hu5", "1hu6", "1hu7", "1f0d", "1f0e", "1f0g")
detailed_data_full <-data.frame()

numberofpdb <- length( pdb_full )

for( index in 1:numberofpdb){
  
  print(index)
  pdb <- pdb_full[index]
  print(pdb)
  
  workdir1 = paste(workdir, "git_test1and2" ,sep="/")
  orientation.frame( df.tm.natives, pdb, workdir1 )
  
  
  df.tm.flag0 <- data.frame( read.table( file = paste( workdir1, paste(pdb,"_min.dat", sep=""), sep = "/"), header=T) )
                     
  bestorientation.df <- compute.best.orientation( df.tm.flag0, 0 )
  min_tilt = bestorientation.df$best.x
  #print(min_tilt)
  min_zcoord = bestorientation.df$best.y
  #print(min_zcoord)
  min_azim = bestorientation.df$best.z
  
  p1 <- plot.orientation.map( workdir1, 0, pdb ) + ggtitle("Minimization on azimuthal angle")
  
  #df.tm.flag1 <- orientation.frame(df.tm.natives, pdb, 1)
  df.tm.flag1 <- data.frame(read.table( file = paste( workdir1, paste(pdb,"_previous.dat", sep=""), sep = "/"), header=T ))
  bestorientation.df <- compute.best.orientation( df.tm.flag1, 1 )
  prev_tilt = bestorientation.df$best.x
  prev_zcoord = bestorientation.df$best.y
  prev_azim = bestorientation.df$best.z
  
  size_diff <- length(prev_tilt) - length(min_tilt)
  max_size <- max( length(prev_tilt), length(min_tilt) )
 
  count_loop <- 0
  
  if( size_diff > 0 )
  {
    while(count_loop < size_diff){
   
      min_tilt <- c( min_tilt,'NULL' )
      min_zcoord <- c( min_zcoord, 'NULL' )
      min_azim <- c( min_azim, 'NULL' )
      count_loop = count_loop+1
     }
    
  } else if( size_diff < 0 )
  {
    while(count_loop < abs(size_diff)){
      
      prev_tilt <- c( prev_tilt,'NULL' )
      prev_zcoord <- c( prev_zcoord, 'NULL' )
      prev_azim <- c( prev_azim, 'NULL' )
      count_loop = count_loop+1
  }
    
    
  }
  
  p2 <- plot.orientation.map( workdir1, 1, pdb) + ggtitle("At azimuthal angle = 0")
  
  exp_tilt = experimental.data[ which( experimental.data$pdb.id == pdb ), 2]
  exp_zcoord = experimental.data[ which( experimental.data$pdb.id == pdb ), 3]
  
    
  detailed_data <- data.frame( min_tilt, min_zcoord, min_azim, prev_tilt, prev_zcoord, prev_azim )
  detailed_data$pdb.id = pdb
  detailed_data$exp_tilt = exp_tilt
  detailed_data$exp_zcoord = exp_zcoord
  
    #print(detailed_data)
  detailed_data_full <- rbind(detailed_data_full, detailed_data)
  #df.tm.flag2 <- orientation.frame(df.tm.natives, pdb, 2)
  p3 <- plot.orientation.map( workdir1, 2, pdb) + ggtitle("At zcoord z = 15")
  #df.tm.flag3 <- orientation.frame(df.tm.natives, pdb, 3)
  p4 <- plot.orientation.map( workdir1, 3, pdb) + ggtitle("At zcoord z = -15")
  
  plot_grid( p1, p2, p3, p4, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 9)
  a1.grid <- plot_grid( p1, p2, p4, p3, ncol = 2, nrow = 2, labels = c("a", "b", "c", "d"), label_size = 9 )
  save_plot( paste( workdir1, paste(pdb,sep="_","full.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )

  plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9)
  a1.grid <- plot_grid( p1, p2, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9 )
  #save_plot( paste( workdir, paste(pdb,sep="_","paperupdates.eps"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )
  save_plot( paste( workdir1, paste(pdb,sep="_","paperupdates.pdf"), sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )
  
  
}
 write.table( detailed_data_full, file = paste( workdir, "processed_data_test1and2.txt", sep = "/") )

