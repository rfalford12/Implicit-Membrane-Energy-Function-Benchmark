# @file: plot_KL_divergence_rates.R
# @brief: Plot KL divergence rates
# @author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(reshape2)

# Current working directory
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/BioInspired-Membrane/figures/figure4-design-results-all"
df <- read.table( paste( workdir, "per_amino_acid_divergence.dat", sep = "/"), header = T )
df$residue <- factor( df$residue, levels = c("A", "I", "L", "M", "V", "F", "W", "Y", "C", "G", "P", "N", "Q", "S", "T", "H", "D", "E", "R", "K"))
df$category <- factor( df$category, levels = c( "nonpolar", "aromatic", "special", "polar", "charged"))

## Just compare overall AA composition statistics
overall.df <- melt( df, ids = c("residue", "category", "efxn"))
p <- ggplot( overall.df[ which( overall.df$variable == "overall" & overall.df$efxn != "m19_all_DLPC" & overall.df$efxn != "r15_all_DLPC" ), ], aes( x = residue, y = value, fill = category ) ) + 
  geom_hline( yintercept = 0, color = "black", size = 0.25 ) + 
  geom_bar( stat = "identity", color = "black", size = 0.15) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  #scale_fill_manual( values = c( "#984ea3", "#e41a1c", "#ffff33", "#377eb8", "#ff7f00") ) + 
  scale_y_continuous( limits = c(-5, 1), breaks = c(-5, -4, -3, -2, -1, 0, 1), "Divergence" ) + 
  scale_x_discrete( "" ) + 
  facet_wrap( ~ efxn, scales = "free" ) + 
  theme_bw() + 
  theme( strip.background = element_blank(), 
         legend.position = "none", 
         text = element_text( size = 8 ), 
         axis.text = element_text( size = 7 ), 
         axis.line.x = element_line( size = 0.25 ), 
         axis.line.y = element_line( size = 0.25 ), 
         axis.ticks = element_line( size = 0.25 ), 
         panel.grid.major = element_line( size = 0.15 ), 
         panel.grid.minor = element_line( size = 0.15 ), 
         panel.border = element_rect( size = 0.25 ) )
print(p)
save_plot( paste( workdir, "divergence_per_AA.png", sep = "/" ), p, units = "in", base_width = 3.42, base_height = 3.42 )

q <- ggplot( overall.df[ which( overall.df$variable == "lipid_facing" & overall.df$efxn != "m19_all_DLPC" ), ], aes( x = residue, y = value, fill = category ) ) + 
  geom_hline( yintercept = 0, color = "black" ) + 
  geom_bar( stat = "identity", color = "black", size = 0.15) + 
  scale_fill_manual( values = c( "#984ea3", "#e41a1c", "#ffff33", "#377eb8", "#ff7f00") ) + 
  scale_y_continuous( limits = c(-16, 5), "Divergence" ) + 
  scale_x_discrete( "" ) + 
  facet_wrap( ~ efxn, scales = "free" ) + 
  theme_bw() + 
  theme( strip.background = element_blank(), 
         legend.position = "none", 
         text = element_text( size = 8 ), 
         axis.text = element_text( size = 7 ), 
         axis.line.x = element_line( size = 0.25 ), 
         axis.line.y = element_line( size = 0.25 ), 
         axis.ticks = element_line( size = 0.25 ), 
         panel.grid.major = element_line( size = 0.25 ), 
         panel.grid.minor = element_line( size = 0.25 ) )
print(q)
save_plot( paste( workdir, "divergence_per_AA.pdf", sep = "/" ), p, units = "in", base_width = 5, base_height = 1.75 )

# Plot the native distribution
df2 <- read.table( paste( workdir, "native_aa_distribution.txt", sep = "/"), header = T) 
df2$Category <- factor( df2$Category, levels = c( "nonpolar", "aromatic", "special", "polar", "charged"))
native.pie <- ggplot( df2, aes( x = "", y = native, fill = Category ) ) + 
  geom_bar( stat = "identity" ) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )

m07.pie <- ggplot( df2, aes( x = "", y = m07, fill = Category ) ) + 
  geom_bar( stat = "identity") + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )

m12.pie <- ggplot( df2, aes( x = "", y = m12, fill = Category ) ) + 
  geom_bar( stat = "identity" ) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank()  )
  
m19.pie <- ggplot( df2, aes( x = "", y = m19, fill = Category ) ) + 
  geom_bar( stat = "identity" ) + 
  scale_fill_brewer( palette = "Pastel1" ) + 
  coord_polar("y", start = 0 ) + 
  theme( legend.position = "none", 
         axis.line = element_blank(), 
         axis.text = element_blank(),
         axis.title = element_blank(), 
         axis.ticks = element_blank() )

combined.pie <- plot_grid( m07.pie, m12.pie, m19.pie, nrow = 1, ncol = 3, scale = c(1.1, 1.1, 1.1) )

save_plot( paste( workdir, "m07_pie.pdf", sep = "/"), m07.pie, units = "in", base_width = 0.5, base_height = 0.5 )
save_plot( paste( workdir, "m12_pie.pdf", sep = "/"), m12.pie, units = "in", base_width = 0.5, base_height = 0.5 )
save_plot( paste( workdir, "m19_pie.pdf", sep = "/"), m19.pie, units = "in", base_width = 0.5, base_height = 0.5 )
# Updated Figure 4
#figure4 <- plot_grid( p, combined.pie, ncol = 1, nrow = 2 )
#save_plot( paste( workdir, "figure4_sample.pdf", sep = "/"), figure4, units = "in", base_width = 5.75, base_height = 4 )

