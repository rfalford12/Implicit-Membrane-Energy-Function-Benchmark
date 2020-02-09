# @file: plot_decoy_discrimination.R
# @author: Rebecca F. Alford (ralford3@jhu.edu)
# @brief: Compare low resolution and high resolution decoy discrimination

library(ggplot2)
library(cowplot)
library(viridis)

workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/decoy-discrimination/data"

read.models <- function( workdir ) {
  df <- data.frame()
  data.files <- list.files(path = workdir, pattern = "*.sc")
  for(f in 1:length(data.files)) {
    full.path <- paste( workdir, data.files[f], sep = "/" )
    sub.df <- read.table( full.path, header = T )
    sub.df$pdb <- unlist(strsplit( data.files[f], "_" )[1])[1]
    sub.df$resolution <- unlist(strsplit( data.files[f], "_" )[1])[2]
    sub.df <- subset( sub.df, select = c("total_score", "rms", "pdb", "resolution" ) )
    df <- rbind( df, sub.df )
  }
  return(df)
}

df <- read.models(workdir)


refined.natives <- read.table( paste( workdir, "refined_natives.txt", sep = "/"), header = T )
combined.df <- rbind( df, refined.natives )
combined.df$resolution <- factor(combined.df$resolution, levels = c("alowres", "hires", "native"))

p <- ggplot( data = combined.df, aes( x = rms, y = total_score, color = resolution ) ) +
  theme_bw() + 
  background_grid() + 
  geom_point( size = 0.35 ) +
  scale_x_continuous( "RMS (Å)", limits = c(0, 40) ) + 
  scale_y_continuous( "Score (REU)", limits = c(NA, -200 ) ) + 
  scale_color_manual( values = c("#377eb8", "#e41a1c", "#4daf4a") ) + 
  facet_wrap( ~ pdb, ncol = 5, nrow = 1, scales = "free_y") + 
  theme( legend.position = "none", 
         text = element_text( size = 8, color = "black" ), 
         axis.text = element_text( size = 8, color = "black" ), 
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35 ), 
         panel.border = element_rect( color = "black" ))

save_plot( paste( workdir, "../", "decoy_disc.pdf", sep = "/"), p, units = "in", base_width = 6.8, base_height = 1.75 )

