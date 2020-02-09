#@file: plot_docked_models.R
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Plot docking models relative to local refined natives

library(cowplot)

# This is temporary
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/protein-protein-docking/D3_heterodimers_bound"

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
  df$CAPRI_rank <- as.factor( df$CAPRI_rank )
  return(df)
  
}

df <- read.homodimer.cases(workdir)

p <- ggplot( data = df, aes( x = Irms, y = I_sc, color = CAPRI_rank ) ) + 
  background_grid() + 
  geom_point() + 
  scale_x_continuous( "RMS (Å)", limits = c(0, 20), expand = c(0,0) ) + 
  scale_y_continuous( "ddG of binding", limits = c(NA, 0) ) + 
  scale_color_manual( values = c("#bdbdbd", "#377eb8", "#e41a1c", "#ff7f00", "#984ea3" ) ) + 
  facet_wrap( ~ pdb.id, scales = "free_y" )
print(p)

