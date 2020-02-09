#@file: plot_docked_models.R
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Plot docking models relative to local refined natives

library(cowplot)

# This is temporary
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/protein-protein-docking/D4_heterodimers_unbound"

read.homodimer.cases <- function(dir) {
  
  df <- data.frame()
  
  docking.model.files <- list.files( dir, pattern = "*_data.sc" )
  for (i in docking.model.files) {
    # Read in docking models
    sub.df <- read.table( paste( dir, i, sep = "/"), header = T, fill = T )
    pdb.id <- unlist(strsplit(i, "_")[1])[1]
    partners <- unlist(strsplit(i, "_")[1])[2]
    sub.df$pdb.id <- pdb.id
    sub.df <- sub.df[ which( sub.df$CAPRI_rank == 0 | sub.df$CAPRI_rank == 1 | sub.df$CAPRI_rank == 2 | sub.df$CAPRI_rank == 3 ), ]
    
    # Read in local models file
    local.models.file = paste( pdb.id, partners, "data", "local", sep = "_")
    local.models.file = paste( local.models.file, "sc", sep = ".")
    local.models <- read.table( paste(dir, local.models.file, sep = "/" ), header = T)
    local.models$pdb.id <- pdb.id
    local.models$CAPRI_rank <- "local"

    # Bind them together
    sub.df <- rbind( sub.df, local.models )
    df <- rbind( df, sub.df )
  }
  #df$CAPRI_rank <- as.factor( df$CAPRI_rank )
  return(df)
  
}

df <- read.homodimer.cases(workdir)
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
save_plot( paste( workdir, "D4_unbound_targets.pdf", sep = "/"), p, units = "in", base_width = 3.5, base_height = 3.3 )

