#@file: plot_docked_models.R
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Plot docking models relative to local refined natives

library(cowplot)
library(boot)

#workdir <- "/Volumes/ralford/Implicit-Membrane-Energy-Function-Benchmark/data/protein-protein-docking/franklin2019/D2_single_TM_complexes/scorefiles"
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/protein-protein-docking/D3_heterodimers_bound"

# Top scoring near-native model count (number of near native decoys
# among a set number of top scoring decoys)
compute.n5 <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc),]
  top.scoring <- head(sorted.df, 5)
  n.near.native <- length(which(top.scoring$CAPRI_rank != 0 ))
  return(n.near.native)
}

compute.n50 <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc),]
  top.scoring <- head(sorted.df, 50)
  n.near.native <- length(which(top.scoring$CAPRI_rank != 0 ))
  return(n.near.native)
}

compute.n100 <- function(df, indices) {
  sub.df <- df[indices,]
  sorted.df <- sub.df[ order(sub.df$I_sc),]
  top.scoring <- head(sorted.df, 100)
  n.near.native <- length(which(top.scoring$CAPRI_rank != 0 ))
  return(n.near.native)
}

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

    # Bind them together
    df <- rbind( df, metric.df )
  }
  return(df)
}

# read.homodimer.cases <- function(dir) {
#   
#   df <- data.frame()
#   
#   docking.model.files <- list.files( dir, pattern = "*_data.sc" )
#   for (i in docking.model.files) {
#     # Read in docking models
#     sub.df <- read.table( paste( dir, i, sep = "/"), header = T )
#     pdb.id <- unlist(strsplit(i, "_")[1])[1]
#     sub.df$pdb.id <- pdb.id
#     
#     # Read in local models file
#     local.models.file = paste( pdb.id, "AB", "data", "local", sep = "_")
#     local.models.file = paste( local.models.file, "sc", sep = ".")
#     local.models <- read.table( paste(dir, local.models.file, sep = "/" ), header = T)
#     local.models$pdb.id <- pdb.id
#     local.models$CAPRI_rank <- "local"
#     
#     # Bind them together
#     sub.df <- rbind( sub.df, local.models )
#     df <- rbind( df, sub.df )
#   }
#   df$CAPRI_rank <- as.factor( df$CAPRI_rank )
#   df <- df[ which( df$CAPRI_rank == 0 | df$CAPRI_rank == 1 | df$CAPRI_rank == 2 | df$CAPRI_rank == 3 ), ]
#   return(df)
#   
# }
# 
# df <- read.homodimer.cases(workdir)

# Compute PNear value for funnels
# For an interface score, the maximum would be zero
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

# Compute docking metrics
# metrics.df <- compute.docking.metrics(workdir)
# pnear.df <- compute.pnear(workdir)
# metrics.df$sampled.rms <- pnear.df$sampled.rms
# metrics.df$pnear <- pnear.df$pnear
# metrics.df$weighted.rms <- pnear.df$weighted.rms
# View(metrics.df)

# p <- ggplot( data = df, aes( x = Irms, y = I_sc, color = CAPRI_rank ) ) + 
#   background_grid() + 
#   geom_point() + 
#   scale_x_continuous( "RMS (Å)", limits = c(0, 30), expand = c(0,0) ) + 
#   scale_y_continuous( "ddG of binding", limits = c(-40, 0) ) + 
#   scale_color_manual( values = c("#bdbdbd", "#377eb8", "#e41a1c", "#ff7f00", "#984ea3" ) ) + 
#   facet_wrap( ~ pdb.id )
# print(p)

# Make the first panel of the docking figure
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/protein-protein-docking"
homodimer.df <- read.table( paste( workdir, "D2_homodimer_metrics.dat", sep = "/"), header = T )
heterodimer.df <- read.table( paste( workdir, "D3_heterodimer_metrics.dat", sep = "/"), header = T )

p <- ggplot() +
  theme_bw() +
  background_grid() +
  geom_vline( xintercept = 1.0, color = "#e41a1c", size = 0.35 ) +
  stat_ecdf( data = homodimer.df, aes( x = e1 ), geom = "density", position = "identity", color = "#737373", size = 0.35 ) +
  stat_ecdf( data = homodimer.df, aes( x = e10 ), geom = "density", position = "identity", color = "#737373", linetype = "longdash", size = 0.35 ) +
  stat_ecdf( data = heterodimer.df, aes( x = e1 ), geom = "density", position = "identity", color = "#000000", size = 0.35 ) +
  stat_ecdf( data = heterodimer.df, aes( x = e10 ), geom = "density", position = "identity", color = "#000000", linetype = "longdash", size = 0.35 ) +
  scale_x_continuous( "<Enrichment>", limits = c(0, NA), expand = c(0,0) ) +
  scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) +
  theme( legend.position = "none",
         text = element_text( size = 8, color = "black" ),
         axis.text = element_text( size = 8, color = "black" ),
         axis.line = element_blank(),
         axis.ticks = element_line( size = 0.35, color = "black" ) )

# Make a merged dataset for histogram
homodimer.df$type <- "homodimers"
heterodimer.df$type <- "heterodimers"
all.df <- rbind(homodimer.df, heterodimer.df)

q <- ggplot() +
  theme_bw() +
  background_grid() +
  geom_vline( xintercept = 3.0, color = "#e41a1c", size = 0.35 ) +
  geom_density( data = all.df, aes( x = n5, color = type ), position = "dodge", size = 0.35 ) +  
  scale_color_manual( values = c("#737373", "#000000") ) +
  scale_x_continuous( "<N5>", limits = c(0, 5), expand = c(0,0) ) +
  scale_y_continuous( "Density", limits = c(0,1.0), expand = c(0,0) ) +
  theme( legend.position = "none",
         text = element_text( size = 8, color = "black" ),
         axis.text = element_text( size = 8, color = "black" ),
         axis.line = element_blank(),
         axis.ticks = element_line( size = 0.35, color = "black" ) )

combined <- plot_grid(p,q, ncol = 1, nrow = 2, labels = c("a", "b"), label_size = 10 )

save_plot( paste( workdir, "D2_and_D3_enrichment.pdf", sep = "/"), combined, units = "in", base_width = 1.85, base_height = 3.3 )



