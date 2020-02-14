# @file: plot_side_chain_distribution.R
# @author: Rebecca F. Alford (ralford3@jhu.edu)
# @brief: Compare desgned and native side chain distribution

library(cowplot)

workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/sc-distribution"

native.distribution <- read.table( paste( workdir, "native_side_chain_distribution.txt", sep = "/"), header = T)
design.distribution <- read.table( paste( workdir, "design_side_chain_distribution.txt", sep = "/"), header = T)
df <- rbind( native.distribution, design.distribution )

# Calculate the difference in density
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

# Calculate the difference in area
get.auc.diff <- function(df, aa) {
  auc.native <- get.auc(df, aa, "native")
  auc.design <- get.auc(df, aa, "design")
  diff <- abs(auc.native - auc.design)
  return(diff)
}

# Do this for the 20 amino acids
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

q <- ggplot( data = df.auc.diff, aes( x = reorder( AA, diffs), y = diffs ) ) + 
  geom_bar( stat = "identity", position = "dodge" ) + 
  scale_x_discrete( "Amino Acids", expand = c(0,0) ) + 
  scale_y_continuous( "AUC", expand = c(0,0), limits = c(0, 0.08) ) + 
  theme_bw() +
  theme( legend.position = "none",
         text = element_text( size = 8, color = "black"),
         axis.text = element_text( size = 8, color = "black" ) )
save_plot( "~/Desktop/auc_diffs.pdf", q, units = "in", base_width = 4, base_height = 3 )
View(df.auc.diff)

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
save_plot( paste( workdir, "figure6_side_chain_distribution.pdf", sep = "/" ), p, units = "in", base_width = 6.85, base_height = 4)


