#@file: plot_helix_kink_data.R
#@author: Rebecca F. Alford (ralford3@jhu.edu)
#@brief: Plot information about kinked helices

library(cowplot)
library(reshape2)
library(ggplot2)

workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/helix-kinks"

# Read in comparison of 2rh1 & 3p0g
gpcr.inactive <- read.table( paste( workdir, "2rh1/output_kinks.txt", sep = "/"), header = T )
gpcr.inactive.scores <- read.table( paste( workdir, "2rh1/2rh1_f19.sc", sep = "/"), header = T )
gpcr.inactive$score <- gpcr.inactive.scores$total_score 
gpcr.inactive.df <- melt( gpcr.inactive, id=c("model", "score") )
gpcr.inactive.df$state <- "inactive"

gpcr.active <- read.table( paste( workdir, "3p0g/output_kinks.txt", sep = "/"), header = T )
gpcr.active.scores <- read.table( paste( workdir, "3p0g/3p0g_f19.sc", sep = "/"), header = T )
gpcr.active$score <- gpcr.active.scores$total_score 
gpcr.active.df <- melt( gpcr.active, id=c("model", "score") )
gpcr.active.df$state <- "active"
gpcr <- rbind( gpcr.active.df, gpcr.inactive.df )

# Read in the native kink angles
native.inactive <- read.table( paste( workdir, "2rh1/output_kinks_native.txt", sep = "/" ), header = T)
native.inactive$state <- "inactive"
native.active <- read.table( paste( workdir, "3p0g/output_kinks_native.txt", sep = "/"), header = T )
native.active$state <- "active"
natives.gpcr <- rbind(native.inactive, native.active)

# Data for potassium receptor (open/closed states)
kcsa.closed <- read.table( paste( workdir, "1r3j/output_kinks.txt", sep = "/"), header = T )
kcsa.closed.scores <- read.table( paste( workdir, "1r3j/1r3j_f19.sc", sep = "/"), header = T )
kcsa.closed$score <- kcsa.closed.scores$total_score 
kcsa.closed.df <- melt(kcsa.closed, id=c("model", "score") )
kcsa.closed.df$state <- "closed"

kcsa.open <- read.table( paste( workdir, "3f7v/output_kinks.txt", sep = "/"), header = T )
kcsa.open.scores <- read.table( paste( workdir, "3f7v/3f7v_f19.sc", sep = "/"), header = T )
kcsa.open$score <- kcsa.open.scores$total_score 
kcsa.open.df <- melt(kcsa.open, id=c("model", "score") )
kcsa.open.df$state <- "open"
kcsa <- rbind( kcsa.closed.df, kcsa.open.df )

# Read in the native kink angles
native.closed <- read.table( paste( workdir, "1r3j/output_kinks_native.txt", sep = "/" ), header = T)
native.closed$state <- "closed"
native.open <- read.table( paste( workdir, "3f7v/output_kinks_native.txt", sep = "/"), header = T )
native.open$state <- "open"
natives.kcsa <- rbind(native.open, native.closed )

# Read in comparison of Adiponectin receptor 1
adipo.closed <- read.table( paste( workdir, "3wxv/output_kinks.txt", sep = "/"), header = T )
adipo.closed.scores <- read.table( paste( workdir, "3wxv/3wxv_f19.sc", sep = "/"), header = T )
adipo.closed$score <- adipo.closed.scores$total_score 
adipo.closed.df <- melt( adipo.closed, id=c("model", "score") )
adipo.closed.df$state <- "closed"

adipo.open <- read.table( paste( workdir, "5lxg/output_kinks.txt", sep = "/"), header = T )
adipo.open.scores <- read.table( paste( workdir, "5lxg/5lxg_f19.sc", sep = "/"), header = T )
adipo.open$score <- adipo.open.scores$total_score 
adipo.open.df <- melt( adipo.open, id=c("model", "score") )
adipo.open.df$state <- "open"
adipo <- rbind( adipo.open.df, adipo.closed.df )

native.closed <- read.table( paste( workdir, "3wxv/output_kinks_native.txt", sep = "/" ), header = T)
native.closed$state <- "closed"
native.open <- read.table( paste( workdir, "5lxg/output_kinks_native.txt", sep = "/"), header = T )
native.open$state <- "open"
natives.adipo <- rbind(native.open, native.closed )

# Read in comparison of Platelet activating receptor
platelet.closed <- read.table( paste( workdir, "5zkp/output_kinks.txt", sep = "/"), header = T )
platelet.closed.scores <- read.table( paste( workdir, "5zkp/5zkp_f19.sc", sep = "/"), header = T )
platelet.closed$score <- platelet.closed.scores$total_score 
platelet.closed.df <- melt( platelet.closed, id=c("model", "score") )
platelet.closed.df$state <- "closed"

platelet.open <- read.table( paste( workdir, "5zkq/output_kinks.txt", sep = "/"), header = T )
platelet.open.scores <- read.table( paste( workdir, "5zkq/5zkq_f19.sc", sep = "/"), header = T )
platelet.open$score <- platelet.open.scores$total_score 
platelet.open.df <- melt( platelet.open, id=c("model", "score") )
platelet.open.df$state <- "open"
platelet <- rbind( platelet.open.df, platelet.closed.df )

native.closed <- read.table( paste( workdir, "5zkp/output_kinks_native.txt", sep = "/" ), header = T)
native.closed$state <- "closed"
native.open <- read.table( paste( workdir, "5zkq/output_kinks_native.txt", sep = "/"), header = T )
native.open$state <- "open"
natives.platelet <- rbind(native.open, native.closed )

p <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  geom_point( data = kcsa[ which( kcsa$variable == "helix3"),], aes( x = value, y = score ), color = "#e41a1c", size = 0.35 ) +
  geom_point( data = kcsa[ which( kcsa$variable == "helix4"),], aes( x = value, y = score ), color = "#253494", size = 0.35 ) +
  geom_vline( data = natives.kcsa, aes( xintercept = helix3 ), color = "#e41a1c", linetype = "dashed") + 
  geom_vline( data = natives.kcsa, aes( xintercept = helix4 ), color = "#253494", linetype = "dashed") + 
  scale_x_continuous( "Kink Angle (degrees)", limits = c(25, 70), expand = c(0,0) ) + 
  scale_y_continuous("Score (REU)" ) + 
  facet_wrap( ~ state, scales = "free_y", ncol = 1, nrow = 2 ) + 
  theme( text = element_text( size = 8 ), 
         axis.text = element_text( size = 8 ), 
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35 )) 
save_plot( paste( workdir, "kcsa_active_v_inactive.pdf", sep = "/"), p, units = "in", base_height = 2.75, base_width = 2 )

q <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  #geom_point( data = gpcr[ which( gpcr$variable == "helix5"),], aes( x = value, y = score ), color = "#e41a1c", size = 0.35) +
  geom_point( data = gpcr[ which( gpcr$variable == "helix6"),], aes( x = value, y = score, color = state ), size = 1.5 ) +
  #geom_point( data = gpcr[ which( gpcr$variable == "helix7"),], aes( x = value, y = score ), color = "#984ea3", size = 0.35 ) +
  #geom_vline( data = natives.gpcr, aes( xintercept = helix5 ), color = "#e41a1c", linetype = "dashed") + 
  geom_vline( data = natives.gpcr, aes( xintercept = helix6, color = state ), linetype = "dashed", size = 1.5 ) + 
  #geom_vline( data = natives.gpcr, aes( xintercept = helix7 ), color = "#984ea3", linetype = "dashed") + 
  scale_color_manual( values = c("#e41a1c", "#253494")) + 
  scale_x_continuous( "Kink Angle (degrees)", limits = c(0, 90), expand = c(0,0) ) + 
  scale_y_continuous("Score (REU)" ) + 
  facet_wrap( ~ state, scales = "free_y", ncol = 2, nrow = 1 ) + 
  theme( legend.position = "none",  
         text = element_text( size = 20, color = "black" ), 
         axis.text = element_text( size = 20, color = "black" ), 
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 1.5), 
         panel.border = element_rect( size = 1.5 )) 
print(q)
save_plot( paste( "~/Desktop", "gpcr_open_v_closed.pdf", sep = "/"), q, units = "in", base_height = 4.5, base_width = 6.5 )

r <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  geom_point( data = adipo[ which( adipo$variable == "helix5"),], aes( x = value, y = score ), color = "#e41a1c", size = 0.35) +
  geom_point( data = adipo[ which( adipo$variable == "helix6"),], aes( x = value, y = score ), color = "#253494", size = 0.35 ) +
  geom_point( data = adipo[ which( adipo$variable == "helix7"),], aes( x = value, y = score ), color = "#984ea3", size = 0.35 ) +
  geom_vline( data = natives.adipo, aes( xintercept = helix5 ), color = "#e41a1c", linetype = "dashed") + 
  geom_vline( data = natives.adipo, aes( xintercept = helix6 ), color = "#253494", linetype = "dashed") + 
  geom_vline( data = natives.adipo, aes( xintercept = helix7 ), color = "#984ea3", linetype = "dashed") + 
  scale_x_continuous( "Kink Angle (degrees)", limits = c(0, 90), expand = c(0,0) ) + 
  scale_y_continuous("Score (REU)" ) + 
  facet_wrap( ~ state, scales = "free_y", ncol = 1, nrow = 2 ) + 
  theme( text = element_text( size = 8 ), 
         axis.text = element_text( size = 8 ), 
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35 )) 
#print(r)
#save_plot( paste( workdir, "adipo_open_v_closed.pdf", sep = "/"), r, units = "in", base_height = 2.75, base_width = 2 )

s <- ggplot() + 
  theme_bw() + 
  background_grid() + 
  geom_point( data = platelet[ which( platelet$variable == "helix5"),], aes( x = value, y = score ), color = "#e41a1c", size = 0.35) +
  geom_point( data = platelet[ which( platelet$variable == "helix6"),], aes( x = value, y = score ), color = "#253494", size = 0.35 ) +
  geom_point( data = platelet[ which( platelet$variable == "helix7"),], aes( x = value, y = score ), color = "#984ea3", size = 0.35 ) +
  geom_vline( data = natives.platelet, aes( xintercept = helix5 ), color = "#e41a1c", linetype = "dashed") + 
  geom_vline( data = natives.platelet, aes( xintercept = helix6 ), color = "#253494", linetype = "dashed") + 
  geom_vline( data = natives.platelet, aes( xintercept = helix7 ), color = "#984ea3", linetype = "dashed") + 
  scale_x_continuous( "Kink Angle (degrees)", limits = c(0, 90), expand = c(0,0) ) + 
  scale_y_continuous("Score (REU)" ) + 
  facet_wrap( ~ state, scales = "free_y", ncol = 1, nrow = 2 ) + 
  theme( text = element_text( size = 8 ), 
         axis.text = element_text( size = 8 ), 
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35 )) 
#print(s)
#save_plot( paste( workdir, "platelet_open_v_closed.pdf", sep = "/"), s, units = "in", base_height = 2.75, base_width = 2 )


