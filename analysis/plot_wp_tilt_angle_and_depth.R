# @file: plot_tilt_angle_and_depth.R
# @author: Rebecca F. Alford (ralford3@jhu.edu)
# @brief: Compare experimental and predicted tilt angle and depth

library(cowplot)
library(nnet)
library(ggrepel)

workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/protein-tilt-angle"
# opm.orientation <- read.table( paste( workdir, "opm_orientation.dat", sep = "/"), header = T )
# 
# compute.best.orientation <- function( df ) { 
#   min.index <- which.min( df$total_score )  
#   best.tilt <- df$angle[ min.index ]
#   best.depth <- df$zcoord[ min.index ]
#   
#   # Make tilt "in-phase"
#   if ( best.tilt > 270 ) {
#      best.tilt = best.tilt - 270
#   } else if ( best.tilt > 180 ) {
#     best.tilt = best.tilt - 180
#   } else if ( best.tilt > 90 ) {
#     best.tilt <- best.tilt - 90
#   }
#   
#   if ( best.tilt > 45 ) {
#     best.tilt <- 90 - best.tilt
#   }
#    
#   df <- data.frame( best.tilt = best.tilt, best.depth = best.depth )
#   return(df)
# }
# 
# read.mappings <- function( dir ) {
#   
#   main.df <- data.frame()
#   data.dirs <- list.dirs(path = dir)
#   for(d in 2:length(data.dirs)) {
#     
#     data.file <- list.files(path = data.dirs[d], pattern = "*.dat")
#     pdb.id <- strsplit(data.file, "_")[1]
#     file.name <- paste(data.dirs[d], data.file, sep = "/")
#     
#     df <- read.table( file.name, header = T )
#     pdb.id <- unlist(pdb.id)[1]
#     orientation.df <- compute.best.orientation(df)
#     orientation.df$pdb.id <- pdb.id
#     orientation.df$opm.tilt <- opm.orientation$tilt[ which( opm.orientation$pdb == pdb.id ) ]
#     orientation.df$opm.depth <- opm.orientation$Z[ which( opm.orientation$pdb == pdb.id ) ]
#     main.df <- rbind( main.df, orientation.df )
#   }
#   return(main.df)
# }
# 
# df <- read.mappings( workdir )
# # Compute the residuals
# depth.lm <- lm(opm.depth ~ best.depth, data = df)
# depth.res <- resid(depth.lm)
# df$depth.residual <- depth.res
# 
# angle.lm <- lm(opm.tilt ~ best.tilt, data = df)
# angle.res <- resid(angle.lm)
# df$angle.residual <- angle.res
# write.table( df, paste( workdir, "protein_position.dat", sep = "/"))

df <- read.table( paste( workdir, "protein_position.dat", sep = "/"), header = T )
df2 <- read.table( paste( workdir, "hydrophobic_thk.dat", sep = "/"), header = T )

# Plot correlation between tilt angles
for.ribbon <- seq(-50, 50, by=1)
p <- ggplot() +
  theme_bw() + 
  background_grid() + 
  geom_abline( size = 0.35 ) +
  geom_ribbon( aes(ymin=for.ribbon-10, ymax=for.ribbon+10, x=for.ribbon), alpha = 0.2) +
  geom_point( data = df, aes( x = opm.tilt, y = best.tilt, color = alpha_or_beta ), size = 0.75  ) + 
  coord_cartesian(ylim = c(0, 45), xlim =c(0, 45) ) +
  scale_x_continuous( "OPM (˚)" ) +
  scale_y_continuous( "Predicted (˚)" ) + 
  scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
  theme( legend.position = "none",
         text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         panel.border = element_rect( color = "black"), 
         axis.ticks = element_line( size = 0.35 ) )

# Plot correlation between membrane depth values
q <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  geom_abline( size = 0.35 ) + 
  geom_ribbon( aes(ymin=for.ribbon-5, ymax=for.ribbon+5, x=for.ribbon), alpha = 0.2) +
  geom_point( data = df, aes( x = abs(opm.depth), y = abs(best.depth), color = alpha_or_beta ), size = 0.75  ) + 
  coord_cartesian(ylim = c(0, 10), xlim =c(0, 10) ) +
  scale_x_continuous( "OPM (Å)" ) +
  scale_y_continuous( "Predicted (Å)" ) +
  scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
  theme( legend.position = "none",
         text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         panel.border = element_rect( color = "black"), 
         axis.ticks = element_line( size = 0.35 ) )

# Plot correlation between hydrophobic thickness values
r <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  geom_abline( size = 0.35 ) + 
  geom_ribbon( aes(ymin=for.ribbon-5, ymax=for.ribbon+5, x=for.ribbon), alpha = 0.2) +
  geom_point( data = df2, aes( x = opm.thk, y = best.thk, color = alpha_or_beta ), size = 0.75  ) + 
  coord_cartesian(ylim = c(15, 41), xlim =c(15,41) ) +
  scale_x_continuous( "OPM (Å)" ) +
  scale_y_continuous( "Predicted (Å)" ) +
  scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
  theme( legend.position = "none",
         text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         panel.border = element_rect( color = "black"), 
         axis.ticks = element_line( size = 0.35 ) )

# Residual Plots
s <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  stat_ecdf( data = df, aes( x = abs(depth.residual)), geom = "density", position = "identity", size = 0.35 ) + 
  stat_ecdf( data = df, aes( x = abs(depth.residual), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) +
  scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
  scale_x_continuous( "Depth Residual (Å)", limits = c(0, NA), expand = c(0,0), breaks = c(0, 5, 10) ) + 
  scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) + 
  theme( legend.position = "none",
         text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         panel.border = element_rect( color = "black"), 
         axis.ticks = element_line( size = 0.35 ) )

t <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  stat_ecdf( data = df, aes( x = abs(angle.residual)), geom = "density", position = "identity", size = 0.35 ) + 
  stat_ecdf( data = df, aes( x = abs(angle.residual), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) +
  scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
  scale_x_continuous( "Angle Residual (˚)", limits = c(0, NA), breaks = c(0, 5, 10, 15, 20, 25), expand = c(0,0) ) + 
  scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) + 
  theme( legend.position = "none",
         text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         panel.border = element_rect( color = "black"), 
         axis.ticks = element_line( size = 0.35 ) )

u <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  stat_ecdf( data = df2, aes( x = abs(thk.residual)), geom = "density", position = "identity", size = 0.35 ) + 
  stat_ecdf( data = df2, aes( x = abs(thk.residual), color = alpha_or_beta ), geom = "density", position = "identity", size = 0.35 ) + 
  scale_color_manual( values = c("#e41a1c", "#377eb8") ) + 
  scale_x_continuous( "Thickness Residual (Å)", limits = c(0, NA), expand = c(0,0) ) + 
  scale_y_continuous( "Frequency", limits = c(0,1.0), expand = c(0,0) ) + 
  theme( legend.position = "none",
         text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         panel.border = element_rect( color = "black"), 
         axis.ticks = element_line( size = 0.35 ) )

# Make large grid for figure
orient.and.thk <- plot_grid(q,p,r,s,t,u, ncol = 3, nrow = 2, labels = c("d", "e", "f", "g", "h", "i"), label_size = 9 )
save_plot( paste( workdir, "orientation_and_thk_plots.pdf", sep = "/"), orient.and.thk, units = "in", base_width = 4.5, base_height = 2.5 )                            

