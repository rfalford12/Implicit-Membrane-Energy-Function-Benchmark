#' Analyze depth-dependent side chain distribution (Test #11)
#' 
#' Computes analysis for depth-dependent side chain distribution including: 
#'    - Depth-dependent density plots for each side chain
#'    - Area under the curve for each side chain
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' and Rituparna Samanta <rsamant2[at]jhu.edu>
#' 

library(cowplot)

#' Read side chain distribution data
#' 
#' @param native.file File with depth-dependent density for native proteins
#' @param design.file File with depth-dependent density for designed proteins
#' 
#' @return Dataframe
#'
#' @example read.sc.distribution.data( "native_side_chain_distribution.txt", "design_side_chain_distribution.txt" )
#' 

read.sc.distribution.data <- function( native.file, design.file ) {
  native.distribution <- read.table( paste( workdir, native.file, sep = "/"), header = T)
  design.distribution <- read.table( paste( workdir, design.file, sep = "/"), header = T)
  df <- rbind( native.distribution, design.distribution )
  return(df)
}

#' Compute the area between the native and design curves using numerical integration methods
#' for natve and design ammino acids
#'
#' @param df Data frame with side chain distribution information
#' @param aa Amino acid of interest
#' 
#' @return AUC value
#'
#' @example get.auc(df, "A")
#' 
get.auc <- function(df, aa) {
  
  d <- density.default(df$zcoord[which(df$AA == aa & df$src == "native")], n = 512, cut = 3 )
  xx <- d$x
  dx <- xx[2L] - xx[1L]
  yy <- d$y
  f.native <- approxfun(xx,yy, yleft = 0, yright = 0)
  
  d.design <- density.default(df$zcoord[which(df$AA == aa & df$src == "design")], n = 512, cut = 3 )
  yy.design <- d.design$y
  xx.design <- d.design$x
  f.design <- approxfun(xx.design,yy.design, yleft=0, yright = 0)
  df.native <-data.frame(xaxis = xx, yaxis = yy, yinterpolated = f.native(xx), ydesigninterpolated = f.design(xx), yerror = yy.error)

  yy.error <- abs(f.native(xx)-f.design(xx))
  f.error <- approxfun(xx,yy.error, yleft = 0, yright = 0)
  
  
  p <- ggplot() + 
    theme_bw() +
    geom_line( data = df.native, aes( x = xaxis, y = yaxis ), size = 0.35 ) + 
    geom_point( data = df.native, aes( x = xaxis, y = yinterpolated ), color = "#2171b5", size = 0.35 )+
    geom_line( data = df.native, aes( x = xaxis, y = ydesigninterpolated ), color = "#9137b8", size = 0.35 )
  
  q <- ggplot() + 
    theme_bw() +
    geom_line( data = df.native, aes( x = xaxis, y = yerror ), size = 0.35 )
    
  
  orient.and.thk <- plot_grid(p,q, ncol = 2, nrow = 1, labels = c("a", "b"), label_size = 9 )
  save_plot( paste( workdir, "density_plots_c.pdf", sep = "/"), orient.and.thk, units = "in", base_width = 4.5, base_height = 2.5 )
  
  
  C <- integrate(f.error, min(xx), max(xx))$value
  print(C)
  
  c.native <- integrate(f.native, min(xx), max(xx))$value
  # c.design <- integrate(f.design, min(xx), max(xx))$value
  print(c.native)
  #this value should be ~1.0 since the histograms are normalized. 
  #  print(c.design)
  #  p.unscaled <- integrate(f, 1, max(xx))$value
  #  p.scaled <- p.unscaled / C
  #  return(p.scaled)
  return(C)
}

#' Compute the difference in area under the native and designed curves
#' 
#' @param df Data frame with side chain distribution information
#' @param aa Amino acid of interest
#' 
#' @return Difference in AUC value
#'
#' @example get.auc.diff(df, "A")
#' 
get.auc.diff <- function(df, aa) {
  auc.native <- get.auc(df, aa, "native")
  auc.design <- get.auc(df, aa, "design")
  diff <- abs(auc.native - auc.design)
  return(diff)
}

#'  Make a data frame with the area between the curves for each amino acid
#' 
#' @param df Data frame with side chain distribution information
#' 
#' @return Data frame with area between the curves
#'
#' @example compute.auc.diff.per.side.chain(df)
#' 
compute.auc.diff.per.side.chain <- function(df) {
  df.auc.diff <- data.frame()
  ala.diff <- get.auc(df, "A")
  cys.diff <- get.auc(df, "C")
  asp.diff <- get.auc(df, "D")
  glu.diff <- get.auc(df, "E")
  phe.diff <- get.auc(df, "F")
  gly.diff <- get.auc(df, "G")
  his.diff <- get.auc(df, "H")
  ile.diff <- get.auc(df, "I")
  lys.diff <- get.auc(df, "K")
  leu.diff <- get.auc(df, "L")
  met.diff <- get.auc(df, "M")
  asn.diff <- get.auc(df, "N")
  pro.diff <- get.auc(df, "P")
  gln.diff <- get.auc(df, "Q")
  arg.diff <- get.auc(df, "R")
  ser.diff <- get.auc(df, "S")
  thr.diff <- get.auc(df, "T")
  val.diff <- get.auc(df, "V")
  trp.diff <- get.auc(df, "W")
  tyr.diff <- get.auc(df, "Y")

  amino.acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  categories <- c( "nonpolar", "special", "charged", "charged", "aromatic", "special", "charged", "nonpolar", "charged", "nonpolar", "nonpolar", "polar", "special", "polar", "charged", "polar", "polar", "nonpolar", "aromatic", "aromatic" )
 
  diffs <- c(ala.diff, cys.diff, asp.diff, glu.diff, phe.diff, gly.diff, his.diff, ile.diff, lys.diff, leu.diff, met.diff, asn.diff, pro.diff, gln.diff, arg.diff, ser.diff, thr.diff, val.diff, trp.diff, tyr.diff )
  
  df.auc.diff <- data.frame( AA = amino.acids, diff = diffs, category=categories )
  
  return(df.auc.diff)
}

#'  Make a plot that ranks the differences in AUC for each amino acid
#' 
#' @param df Data frame with AUC differences; area between the curves
#' 
#' @return Plot object
#'
#' @example make.auc.diff.plot(df)
#' 
make.auc.diff.plot <- function(df) {
  q <- ggplot( data = df, aes( x = reorder( AA, diff), y = diff, fill=category ) ) + 
    geom_bar( stat = "identity", position = "dodge" ) + 
    scale_fill_brewer( palette = "Pastel1" ) +
    scale_x_discrete( "Amino Acids", expand = c(0,0) ) + 
    scale_y_continuous( "AUC", expand = c(0,0), limits = c(0, 0.30) ) + 
    theme_bw() +
    theme( legend.position = "top",
           text = element_text( size = 8, color = "black"),
           axis.text = element_text( size = 8, color = "black" ) ) 
  
  return(q)
}

#' Make a plot with the depth-dependent side chain density for each
#' amino acid of the 20 canonical amino acids 
#' 
#' @param df Data frame side chain distributions 
#' @param df1 Data frame with AUC differences; area between the curves
#'
#' @return Plot object
#'
#' @example make.side.chain.density.plot(df, df1)
#' 
make.side.chain.density.plot <- function(df, df1) {
  dataset <- data.frame()
  
  for( f in 1:nrow(df1) ){
    print(as.character(df1[["AA"]])[f])
    subdf <- df[which(df$AA == as.character(df1[["AA"]])[f]), ]
    subdf$diff <-df1[["diff"]][f]
    dataset <- rbind(dataset, subdf)
    
    }
  
   p <- ggplot( data = dataset, aes( x = zcoord, linetype = src ) ) +
    background_grid() +
    stat_density(position="identity", geom="line") +
    scale_color_manual( values = c("#e41a1c", "#377eb8") ) +
    scale_x_continuous( "Membrane Depth (Å)", limits = c(-30, 30), breaks = c(-20, -10, 0, 10, 20), expand = c(0,0) ) +
    scale_y_continuous( "Density" ) +
    facet_wrap( ~ AA, ncol=4 ) +
  #  geom_text( size = 3, mapping = aes(x = -30, y = -Inf, label = round(diff,3)), hjust = -2.5, vjust = -7.5, fontface="plain")+
    theme_bw() +
    theme( legend.position = "none",
           text = element_text( size = 10, color = "black"),
           axis.text = element_text( size = 10, color = "black" ) )
  
  #p<-p+geom_text( data = df1, mapping = aes(x = -Inf, y = -Inf, label = diff), hjust = -0.1, vjust = -1 ) 
  return(p)
}

workdir <- "/path-to-/sequence-recovery"
native.file <- "native_side_chain_distribution.txt"
design.file <- "design_side_chain_distribution.txt"

main.df <- read.sc.distribution.data( native.file, design.file )
diff.df <- compute.auc.diff.per.side.chain( main.df )
write.table( diff.df, file = paste( workdir, "aucdiff_data.txt", sep = "/") )
p1 <- make.auc.diff.plot( diff.df)
a1.grid <- plot_grid( p1, ncol = 1, nrow = 1, label_size = 9 )
save_plot( paste( workdir, "amminoacid_aucdiff_figure_v2.pdf", sep = "/"), a1.grid, units = "in", base_width = 4, base_height = 4 )
p2 <- make.side.chain.density.plot( main.df, diff.df )
a1.grid <- plot_grid( p2, ncol = 1, nrow = 1, label_size = 9 )
save_plot( paste( workdir, "amminoacid_depthwise_figure_v2.pdf", sep = "/"), a1.grid, units = "in", base_width = 9, base_height = 6 )

