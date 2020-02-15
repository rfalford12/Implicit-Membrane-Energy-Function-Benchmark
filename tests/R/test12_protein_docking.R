#' Analyze membrane protein-protein docking data (Test #12)
#' 
#' Computes analysis for protein-protein docking calculations, such as: 
#'    - Enrichment of acceptable models in the top 1% and 10% scoring models
#'    - Score vs. RMSD analysis, coded by CAPRI rank
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

library(cowplot)
library(boot)

#' Number of acceptable models in the top five ranking models
#'
#' @param df Dataframe containing columns CAPRI_rank
#' @param indices Number of elements to sample for bootstrapping
#'
#' @return N5
#'
#' @examples
#' compute.n5(df, R=1000)
#' 
compute.n5 <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc),]
  top.scoring <- head(sorted.df, 5)
  n.near.native <- length(which(top.scoring$CAPRI_rank != 0 ))
  return(n.near.native)
}

#' Number of acceptable models in the top fifty ranking models
#'
#' @param df Dataframe containing columns CAPRI_rank
#' @param indices Number of elements to sample for bootstrapping
#'
#' @return N50
#'
#' @examples
#' compute.n50(df, R=1000)
#' 
compute.n50 <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc),]
  top.scoring <- head(sorted.df, 50)
  n.near.native <- length(which(top.scoring$CAPRI_rank != 0 ))
  return(n.near.native)
}

#' Number of acceptable models in the top 100 ranking models
#'
#' @param df Dataframe containing columns CAPRI_rank
#' @param indices Number of elements to sample for bootstrapping
#'
#' @return N100
#'
#' @examples
#' compute.n100(df, R=1000)
#' 
compute.n100 <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc),]
  top.scoring <- head(sorted.df, 100)
  n.near.native <- length(which(top.scoring$CAPRI_rank != 0 ))
  return(n.near.native)
}

#' Enrichment for # of acceptable models in the top 1% scoring models
#'
#' @param df Dataframe containing columns CAPRI_rank
#' @param indices Number of elements to sample for bootstrapping
#'
#' @return E_1%
#'
#' @examples
#' compute.1percent.enrich(df, R=1000)
#' 
compute.1percent.enrich <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc), ]
  n.decoys <- nrow(sorted.df)
  n.near.native <- nrow(sorted.df[ which(sorted.df$CAPRI_rank != 0 ), ] )
  n.decoys.in.top.n.percent <- round( 0.01 * nrow(sorted.df) )
  top.n.percent <- head(sorted.df, n.decoys.in.top.n.percent)
  near.native.in.top.n.percent <- length(which(top.n.percent$CAPRI_rank != 0 ))
  numerator <- near.native.in.top.n.percent/n.decoys.in.top.n.percent
  denominator <- n.near.native/n.decoys
  return(numerator/denominator)
}

#' Enrichment for # of acceptable models in the top 10% scoring models
#'
#' @param df Dataframe containing columns CAPRI_rank
#' @param indices Number of elements to sample for bootstrapping
#'
#' @return E_10%
#'
#' @examples
#' compute.10percent.enrich(df, R=1000)
#' 
compute.10percent.enrich <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc), ]
  n.decoys <- nrow(sorted.df)
  n.near.native <- nrow(sorted.df[ which(sorted.df$CAPRI_rank != 0 ), ] )
  n.decoys.in.top.n.percent <- round( 0.1 * nrow(sorted.df) )
  top.n.percent <- head(sorted.df, n.decoys.in.top.n.percent)
  near.native.in.top.n.percent <- length(which(top.n.percent$CAPRI_rank != 0 ))
  numerator <- near.native.in.top.n.percent/n.decoys.in.top.n.percent
  denominator <- n.near.native/n.decoys
  return(numerator/denominator)
}

#' Compute PNear - the likelihood that a model will be near-native
#'
#' @param workdir Working directory with docking score files 
#'
#' @return Dataframe of PNear metrics per target
#'
#' @examples
#' compute.docking.metrics(workdir)
#' 
compute.pnear <- function(dir) {
  
  df <- data.frame()
  
  docking.model.files <- list.files( dir, pattern = "*_data.sc" )
  for (i in docking.model.files) {
    
    # Read in docking models
    sub.df <- read.table( paste( dir, i, sep = "/"), header = T, fill = T )
    pdb.id <- unlist(strsplit(i, "_")[1])[1]
    sub.df$pdb.id <- pdb.id
    print(pdb.id)
    
    sorted.df <- sub.df[ order(sub.df$Irms), ]
    top_2percent <- 0.02 * nrow(sorted.df)
    num_sum = 0.0
    denom_sum = 0.0
    weighted_rms = 0.0
    for (r in 1:nrow(sorted.df)) {
      I_sc <- sorted.df$I_sc[r]
      Irms <- sorted.df$Irms[r]
      energy <- exp(-(I_sc)/0.5)
      weighted_energy <- exp( -(Irms*Irms)/(2*2))*energy
      weighted_rms = weighted_rms + (Irms*energy)
      num_sum = num_sum + weighted_energy
      denom_sum = denom_sum + energy
    }
    
    # Print final values
    sampledRMS <- sorted.df$Irms[tail(top_2percent, n=1)]
    pnear <- num_sum/denom_sum
    weightedRMS <- weighted_rms/denom_sum
    pnear.df <- data.frame( pdb.id = pdb.id, sampled.rms <- sampledRMS, pnear = pnear, weighted.rms = weightedRMS )
    df <- rbind( df, pnear.df )
  }
  return(df)
}


#' Compute all enrichment and top model metrics for a dataset
#'
#' @param workdir Working directory with docking score files 
#'
#' @return Dataframe of metrics 
#'
#' @examples
#' compute.docking.metrics(workdir)
#' 
compute.docking.metrics <- function(dir) {

  df <- data.frame()

  docking.model.files <- list.files( dir, pattern = "*_data.sc" )
  for (i in docking.model.files) {

    # Read in docking models
    sub.df <- read.table( paste( dir, i, sep = "/"), header = T, fill = T )
    pdb.id <- unlist(strsplit(i, "_")[1])[1]
    sub.df$pdb.id <- pdb.id

    print(pdb.id)

    # Calculate docking metrics
    e.1.percent <- boot( data = sub.df, statistic=compute.1percent.enrich, R=1000 )
    e.10.percent <- boot( data = sub.df, statistic=compute.10percent.enrich, R=1000)
    n.5 <- boot(data=sub.df, statistic=compute.n5, R=1000 )
    n.50 <-  boot(data=sub.df, statistic=compute.n50, R=1000 )
    n.100 <-  boot(data=sub.df, statistic=compute.n100, R=1000 )

    metric.df <- data.frame( pdb.id = pdb.id,  e.1 = e.1.percent[1], e.10 = e.10.percent[1], n.5 = n.5[1], n.50 = n.50[1], n100 = n.100[1] )
    df <- rbind( df, metric.df )
  }
  return(df)
}

#' Make a dataframe of all I_sc and I_rms values for all docked models and
#' locally refined native models for each target
#'
#' @param workdir Working directory with docking score files 
#'
#' @return Dataframe of scores and rmsd values for each target 
#'
#' @examples
#' read.homodimer.cases(workdir)
#' 
read.homodimer.cases <- function(dir) {

  df <- data.frame()

  docking.model.files <- list.files( dir, pattern = "*_data.sc" )
  for (i in docking.model.files) {
    # Read in docking models
    sub.df <- read.table( paste( dir, i, sep = "/"), header = T )
    pdb.id <- unlist(strsplit(i, "_")[1])[1]
    sub.df$pdb.id <- pdb.id

    # Read in local models file
    local.models.file = paste( pdb.id, "AB", "data", "local", sep = "_")
    local.models.file = paste( local.models.file, "sc", sep = ".")
    local.models <- read.table( paste(dir, local.models.file, sep = "/" ), header = T)
    local.models$pdb.id <- pdb.id
    local.models$CAPRI_rank <- "local"

    # Bind them together
    sub.df <- rbind( sub.df, local.models )
    df <- rbind( df, sub.df )
  }
  df$CAPRI_rank <- as.factor( df$CAPRI_rank )
  df <- df[ which( df$CAPRI_rank == 0 | df$CAPRI_rank == 1 | df$CAPRI_rank == 2 | df$CAPRI_rank == 3 ), ]
  return(df)

}

#' Make a plot of the cummulative distribution of enrichment values
#'
#' @param workdir  Dataframe with docking metrics 
#'
#' @return ggplot2 object
#'
#' @examples
#' make.ecdf.plot(df)
#' 
make.ecdf.plot <- function(df) {
  p <- ggplot() +
    theme_bw() +
    background_grid() +
    geom_vline( xintercept = 1.0, color = "#e41a1c", size = 0.35 ) +
    stat_ecdf( data = df, aes( x = e1 ), geom = "density", position = "identity", color = "#737373", size = 0.35 ) +
    stat_ecdf( data = df, aes( x = e10 ), geom = "density", position = "identity", color = "#737373", linetype = "longdash", size = 0.35 ) +
    scale_x_continuous( "<Enrichment>", limits = c(0, NA), expand = c(0,0) ) +
    scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) +
    theme( legend.position = "none",
           text = element_text( size = 8, color = "black" ),
           axis.text = element_text( size = 8, color = "black" ),
           axis.line = element_blank(),
           axis.ticks = element_line( size = 0.35, color = "black" ) )
  return(p)
}

#' Make a plot of the I_sc vs. I_rms (or dG_separated) colored by CAPRI_rank
#'
#' @param workdir  Dataframe with score and rms values for each target 
#'
#' @return ggplot2 object
#'
#' @examples
#' make.score.vs.rmsd.plot(df)
#' 
make.score.vs.rmsd.plot <- function(df) {
  
  df$CAPRI_rank[ which( df$CAPRI_rank == 0 ) ] <- "incorrect"
  df$CAPRI_rank[ which( df$CAPRI_rank == 1 ) ] <- "acceptable"
  df$CAPRI_rank[ which( df$CAPRI_rank == 2 ) ] <- "medium"
  df$CAPRI_rank[ which( df$CAPRI_rank == 3 ) ] <- "high"
  df$CAPRI_rank[ which( df$CAPRI_rank == "local" ) ] <- "refined"
  df$CAPRI_rank <- factor( df$CAPRI_rank, levels = c("incorrect", "acceptable", "medium", "high", "refined") )
  
  p <- ggplot( data = df, aes( x = Irms, y = I_sc, color = CAPRI_rank ) ) + 
    theme_bw() + 
    background_grid() +
    geom_point( size = 0.1 ) +
    scale_x_continuous( "RMS (Å)", limits = c(0, 10), breaks = c(0, 5, 10), expand = c(0,0) ) +
    scale_y_continuous( "∆∆G of Binding", limits = c(NA, 0) ) +
    scale_color_manual( "CAPRI Rank", values = c("#bdbdbd", "#377eb8", "#ff7f00", "#e41a1c", "#984ea3" ) ) +
    facet_wrap( ~ pdb.id, scales = "free_y" ) + 
    theme( legend.position = "bottom", 
           text = element_text( size = 8, color = "black" ), 
           axis.text = element_text( size = 8, color = "black" ), 
           axis.line = element_blank(), 
           axis.ticks = element_line( color = "black", size = 0.35 ), 
           panel.border = element_rect( color = "black" ), 
           legend.key.size = unit(0.5, "lines"), 
           legend.text = element_text( size = 8, color = "black"),
           legend.margin = margin(0.05, 0.05, 0.05, 0.05, "cm"), 
           legend.box.margin=margin(-7,12,0,-10))
  
  return(p)
}

