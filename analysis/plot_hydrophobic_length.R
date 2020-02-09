# @file: plot_hydrophobic_length.R
# @author: Rebecca F. Alford (ralford3@jhu.edu)
# @brief: Compare experimental, predicted, and calculated hydrophobic thickness values

library(cowplot)
library(nnet)
library(ggrepel)
library(mokken)
library(pracma)


workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/hydrophobic-length"
exp.data <- read.table( paste( workdir, "experimental_info.dat", sep = "/"), header = T )

exp.avail <- c("1hxx", "1fep", "1nqe", "1qj8", "1qfg", "1qd6", "1qjp", "1m01", "1gzm", "1rzh", "1v55", "1r3j", "2oar", "1yce" )
read.hydrophobic.thk <- function( dir ) {
  
  main.df <- data.frame()
  data.dirs <- list.dirs(path = dir)
  for(d in 2:length(data.dirs)) {
    
    data.file <- list.files(path = data.dirs[d], pattern = "*.dat")
    pdb.id <- strsplit(data.file, "_")[1]
    file.name <- paste(data.dirs[d], data.file, sep = "/")

    df <- read.table( file.name, header = T )
    df$pdb.id <- unlist(pdb.id)[1]
    df$max <- df$thickness[ which.is.max(df$score) ]
    df$normalized <- df$score/max(df$score)
    df$deriv <- gradient(df$score, h1 = df$thickness)
    df$exp.avail <- unlist(pdb.id)[1] %in% exp.avail
    df$score.of.best.thk <- max(df$normalized)
    df$best.thk <- df$thickness[ which.is.max(df$normalized) ]
    # exp.thk <- exp.data$exp[ which(exp.data$pdb == unlist(pdb.id)[1]) ]
    # opm.thk <- exp.data$opm[ which(exp.data$pdb == unlist(pdb.id)[1]) ]
    # new.df <- data.frame( pdb.id = unlist(pdb.id)[1], best.thk = best.thk, exp.thk = exp.thk, opm.thk = opm.thk )
    # 
    main.df <- rbind( main.df, df )
  }
  return(main.df)
}

read.pred.v.experiment <- function( dir ) {
  
  main.df <- data.frame()
  data.dirs <- list.dirs(path = dir)
  for(d in 2:length(data.dirs)) {
    
    data.file <- list.files(path = data.dirs[d], pattern = "*.dat")
    pdb.id <- strsplit(data.file, "_")[1]
    file.name <- paste(data.dirs[d], data.file, sep = "/")
    
    df <- read.table( file.name, header = T )
    df$max <- df$thickness[ which.is.max(df$score) ]
    df$normalized <- df$score/max(df$score)
    
    add.row <- data.frame( pdb.id = unlist(pdb.id)[1], 
                           exp.avail = unlist(pdb.id)[1] %in% exp.avail, 
                           score.of.best.thik = max(df$normalized), 
                           best.thk = df$thickness[ which.is.max(df$normalized) ], 
                           exp.thk = exp.data$exp[ which(exp.data$pdb == unlist(pdb.id)[1]) ], 
                           opm.thk = exp.data$opm[ which(exp.data$pdb == unlist(pdb.id)[1]) ] )
    
    main.df <- rbind( main.df, add.row )
  }
  return(main.df)
}

df <- read.hydrophobic.thk(workdir)
# thk.lm <- lm(opm.thk ~ best.thk, data = df)
# thk.res <- resid(thk.lm)
# df$thk.residual <- thk.res
# write.table( df, paste( workdir, "hydrophobic_thk.dat", sep = "/"))


# Data summary
p <- ggplot( data = df  ) +
  theme_bw() + 
  background_grid() +
  geom_vline( xintercept = 15 ) +
  geom_line( aes( x = thickness, y = normalized, color = pdb.id ) ) +
  geom_point( aes( x = best.thk, y = score.of.best.thk, color = pdb.id ) ) +
  scale_x_continuous( "Thickness (Å)", limits = c(0, 40), expand = c(0,0) ) + 
  scale_y_continuous( "Normalized Score", limits = c(1.0, 1.4), expand  = c(0,0) ) + 
  scale_color_discrete( "PDB" ) + 
  theme( text = element_text( size = 8 ), 
         axis.text = element_text( size = 8 ) )
print(p)
save_plot( paste( workdir, "hydrophobic_legth_trace3.pdf", sep = "/"), p, units = "in", base_height = 4.75, base_width = 6 )

# df2 <- read.pred.v.experiment(workdir)
# 
# 
# # Experimental vs. predicted
# for.ribbon <- seq(0, 50, by=1)
# q <- ggplot() +
#   background_grid() +
#   geom_abline() +
#   geom_ribbon( aes(ymin=for.ribbon-5, ymax=for.ribbon+5, x=for.ribbon), alpha = 0.2) +
#   geom_point( data = df2, aes( x = exp.thk, y = best.thk) ) +
#   #geom_text_repel( data = df2[ which( df2$exp.avail == TRUE ), ], aes( x = exp.thk, y = best.thk, label = pdb.id ) ) +
#   coord_cartesian(ylim = c(15, 41), xlim =c(15,41) ) +
#   scale_x_continuous( "Experiment (Å)" ) +
#   scale_y_continuous( "Predicted (Å)" ) 
# 
# # OPM vs. predicted
# r <- ggplot() +
#   background_grid() +
#   geom_abline() +
#   geom_ribbon( aes(ymin=for.ribbon-5, ymax=for.ribbon+5, x=for.ribbon), alpha = 0.2) +
#   geom_point( data = df2, aes( x = opm.thk, y = best.thk ) ) +
#   geom_text_repel( data = df2, aes( x = opm.thk, y = best.thk, label = pdb.id ) ) +
#   coord_cartesian(ylim = c(15, 41), xlim =c(15,41) ) +
#   scale_x_continuous( "OPM (Å)" ) +
#   scale_y_continuous( "Predicted (Å)" ) 
# 
# # Compute the residuals
# # thk.lm <- lm(opm.thk ~ best.thk, data = df2)
# # thk.res <- resid(thk.lm)
# # df2$residual <- thk.res
# # View(df2)
# # 
# # s <- ggplot( data = df2, aes( x = opm.thk, y = residual, label = pdb.id) ) + 
# #   geom_label() + 
# #   geom_hline( yintercept = -5 ) + 
# #   geom_hline( yintercept = 2.5 ) +
# #   geom_hline( yintercept = -2.5 ) + 
# #   geom_hline( yintercept = 5 )
# # print(s)
# 
# combo.fig <- plot_grid(q,r, ncol = 1, nrow = 2)
# save_plot( paste( workdir, "figure3_hydrophobic_length.pdf", sep = "/"), combo.fig, units = "in", base_width = 3.5, base_height = 6 )
