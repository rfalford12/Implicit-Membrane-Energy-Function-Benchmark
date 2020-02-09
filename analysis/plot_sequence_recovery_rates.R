# @file: plot_sequence_recovery_rates.R
# @brief: Plot sequence recovery rates vs. probability above random recovery
# @author: Rebecca F. Alford (ralford3@jhu.edu)

library(cowplot)
library(ggrepel)

# Current working directory
workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/BioInspired-Membrane/figures/figure4-design-results-all"
df <- read.table( paste( workdir, "design_seqrecov_performance.txt", sep = "/"), header = T )
df$solvation <- factor( df$solvation, c("overall", "buried", "surface", "lipid_facing", "interfacial", "aqueous"))
df$which_panel <- factor( df$which_panel, c("FALSE", "burial", "memb"))
recovery.rates <- ggplot( data = df, aes( x = recovery, y = frac.recov.above.random*20, shape = solvation ) ) +
    background_grid() +
    geom_point(size = 1.5, color = "gray40") + 
    geom_text_repel( aes( label = efxn ), size = 2, force = 2 ) +
    scale_x_continuous( "Recovery (%)", limits = c(0, 0.5), expand = c(0,0) ) +
    scale_y_continuous( "# AAs with non-random Recovery", limits = c(0, 20), expand = c(0,0) ) +
    scale_shape_manual( values = c(15, 16, 1, 2, 17) ) + 
    #scale_fill_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
    #scale_color_manual( values = c( "#b3cde3", "#bdbdbd", "#fbb4ae" ) ) +
    facet_wrap(  ~ which_panel, ncol = 3, nrow = 1, scales = "free" ) +
    theme_bw() + 
    theme( legend.position = "none",
           text = element_text( size = 8 ),
           axis.text = element_text( size = 8 ),
           axis.line.x = element_line( size = 0.25 ),
           axis.line.y = element_line( size = 0.25 ),
           axis.ticks = element_line( size = 0.25 ),
           strip.background = element_blank(), 
           panel.grid.major = element_line( size = 0.15 ), 
           panel.grid.minor = element_line( size = 0.15 ), 
           panel.border = element_rect( size = 0.25 ) ) + 
    coord_flip()

print(recovery.rates)
save_plot( paste( workdir, "recovery_rates.pdf", sep = "/"), recovery.rates, units = "in", base_width = 6.5, base_height = 1.5)
