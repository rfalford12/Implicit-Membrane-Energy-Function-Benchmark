#' Analyze depth-dependent side chain distribution (Test #11)
#' 
#' Computes analysis for depth-dependent side chain distribution including: 
#'    - Depth-dependent density plots for each side chain
#'    - Area under the curve for each side chain
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

library(cowplot)

#' Read side chain distribution data
#' 
#' @param native.file File with depth-dependent density for native proteins
#' @param design.file File with depth-dependent density for designed proteins
#' 
#' @return Dataframe
#'
#' @example read.sc.distribution.data( "natives.txt", "designed.txt" )
#' 
read.sc.distribution.data <- function( native.file, design.file ) {
  native.distribution <- read.table( paste( workdir, native.file, sep = "/"), header = T)
  design.distribution <- read.table( paste( workdir, design.file, sep = "/"), header = T)
  df <- rbind( native.distribution, design.distribution )
  return(df)
}

#' Compute the area under the curve using an integration method
#' 
#' @param df Data frame with side chain distribution information
#' @param aa Amino acid of interest
#' @param src Native or designed? 
#' 
#' @return AUC value
#'
#' @example get.auc(df, "A", "native")
#' 
get.auc <- function(df, aa, src) {
  d <- density.default(df$zcoord[which(df$AA == aa & df$src == src )], n = 512, cut = 3 )
  xx <- d$x
  dx <- xx[2L] - xx[1L]
  yy <- d$y
  f <- approxfun(xx,yy)
  C <- integrate(f, min(xx), max(xx))$value
  p.unscaled <- integrate(f, 1, max(xx))$value
  p.scaled <- p.unscaled / C
  return(p.scaled)
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

#'  Make a data frame with the difference in AUC for each amino acid
#' 
#' @param df Data frame with side chain distribution information
#' 
#' @return Data frame with difference in AUC values
#'
#' @example compute.auc.diff.per.side.chain(df)
#' 
compute.auc.diff.per.side.chain <- function(df) {
  df.auc.diff <- data.frame()
  ala.diff <- get.auc.diff(df, "A")
  cys.diff <- get.auc.diff(df, "C")
  asp.diff <- get.auc.diff(df, "D")
  glu.diff <- get.auc.diff(df, "E")
  phe.diff <- get.auc.diff(df, "F")
  gly.diff <- get.auc.diff(df, "G")
  his.diff <- get.auc.diff(df, "H")
  ile.diff <- get.auc.diff(df, "I")
  lys.diff <- get.auc.diff(df, "K")
  leu.diff <- get.auc.diff(df, "L")
  met.diff <- get.auc.diff(df, "M")
  asn.diff <- get.auc.diff(df, "N")
  pro.diff <- get.auc.diff(df, "P")
  gln.diff <- get.auc.diff(df, "Q")
  arg.diff <- get.auc.diff(df, "R")
  ser.diff <- get.auc.diff(df, "S")
  thr.diff <- get.auc.diff(df, "T")
  val.diff <- get.auc.diff(df, "V")
  trp.diff <- get.auc.diff(df, "W")
  tyr.diff <- get.auc.diff(df, "Y")
  amino.acids <- c("A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y")
  diffs <- c(ala.diff, cys.diff, asp.diff, glu.diff, phe.diff, gly.diff, his.diff, ile.diff, lys.diff, leu.diff, met.diff, asn.diff, pro.diff, gln.diff, arg.diff, ser.diff, thr.diff, val.diff, trp.diff, tyr.diff )
  df.auc.diff <- data.frame( AA = amino.acids, diff = diffs )
}

#'  Make a plot that ranks the differences in AUC for each amino acid
#' 
#' @param df Data frame with AUC differences
#' 
#' @return Plot object
#'
#' @example make.auc.diff.plot(df)
#' 
make.auc.diff.plot <- function(df) {
  q <- ggplot( data = df, aes( x = reorder( AA, diffs), y = diffs ) ) + 
    geom_bar( stat = "identity", position = "dodge" ) + 
    scale_x_discrete( "Amino Acids", expand = c(0,0) ) + 
    scale_y_continuous( "AUC", expand = c(0,0), limits = c(0, 0.08) ) + 
    theme_bw() +
    theme( legend.position = "none",
           text = element_text( size = 8, color = "black"),
           axis.text = element_text( size = 8, color = "black" ) )  
}

#' Make a plot with the depth-dependent side chain density for each
#' amino acid of the 20 canonical amino acids 
#' 
#' @param df Data frame side chain distributions 
#' 
#' @return Plot object
#'
#' @example make.side.chain.density.plot(df)
#' 
make.side.chain.density.plot <- function(df) {
  p <- ggplot( data = df, aes( x = zcoord, linetype = src ) ) +
    background_grid() +
    stat_density(position="identity", geom="line") +
    scale_color_manual( values = c("#e41a1c", "#377eb8") ) +
    scale_x_continuous( "Membrane Depth (Å)", limits = c(-30, 30), breaks = c(-20, -10, 0, 10, 20), expand = c(0,0) ) +
    scale_y_continuous( "Density" ) +
    facet_wrap( ~ AA ) +
    theme_bw() +
    theme( legend.position = "none",
           text = element_text( size = 8, color = "black"),
           axis.text = element_text( size = 8, color = "black" ) )
}
