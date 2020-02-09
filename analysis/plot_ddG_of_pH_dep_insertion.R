# @file: plot_ddG_of_pH_insertion.R
# @author: Rebecca F. Alford (ralford3@jhu.edu)
# @brief: Compare experimental and predicted pH-dependent ddG of insertion values

library(cowplot)
library(viridis)
library(ggrepel)

workdir <- "/Users/ralford/Dropbox/rebecca-research-shared/writing/papers/Membrane-Efxn-Benchmarks/data/ddG-of-pH-insertion/C5_pHLIP_helical_peptides"
exp.data <- read.table( paste( workdir, "pH-inserted-peptides-exp-data.dat", sep = "/"), header = T )

# compute.best.membrane.orientation <- function(df) {
# 
#   # Make a dataframe with only the membrane positions
#   memb.df <- df[ which( abs(df$zcoord) <= 15 ), ]
#   ordered.memb.df <- memb.df[ order(memb.df$total_score), ]
#   head <- head(ordered.memb.df,1)
#   best.memb.tilt <- head$angle
#   best.memb.depth <- head$zcoord
#   best.memb.score <- head$total_score
# 
#   # Create a dataframe that stores all of this information and return it
#   out.df <- data.frame( best.memb.tilt = best.memb.tilt,
#                         best.memb.depth = best.memb.depth,
#                         best.memb.score = best.memb.score )
#   return(out.df)
# }
# 
# compute.best.soluble.orientation <- function(df) {
# 
#   # Make a dataframe with only the membrane positions
#   memb.df <- df[ which( abs(df$zcoord) > 30 ), ]
#   ordered.memb.df <- memb.df[ order(memb.df$total_score), ]
#   head <- head(ordered.memb.df,1)
#   best.soluble.tilt <- head$angle
#   best.soluble.depth <- head$zcoord
#   best.soluble.score <- head$total_score
# 
#   # Create a dataframe that stores all of this information and return it
#   out.df <- data.frame( best.soluble.tilt = best.soluble.tilt,
#                         best.soluble.depth = best.soluble.depth,
#                         best.soluble.score = best.soluble.score )
#   return(out.df)
# }
# 
# read.mappings <- function( dir ) {
# 
#   main.df <- data.frame()
#   data.dirs <- list.dirs(path = dir)
#   for(d in 2:length(data.dirs)) {
# 
#     # Get the file
#     data.file <- list.files(path = data.dirs[d], pattern = "*.dat")
#     file.name <- paste(data.dirs[d], data.file, sep = "/")
# 
#     # Get the file metadata
#     file.info <- strsplit(data.file, "_")[1]
#     pdb.id <- unlist(file.info)[1]
#     pH.value <- unlist(strsplit(data.dirs[d], "_")[1])[5]
# 
#     df <- read.table( file.name, header = T )
# 
#     dG.pH4 <- 0
#     dG.pH8 <- 0
# 
#     if ( pH.value == 4 ) {
#       membrane.orientation <- compute.best.membrane.orientation(df)
#       dG.pH4 <- membrane.orientation$best.memb.score
#       pH4.df <- data.frame( pdb.id = pdb.id, pH.value = pH.value, dG = dG.pH4 )
#       main.df <- rbind( main.df, pH4.df )
#     } else {
#       soluble.orientation <- compute.best.soluble.orientation(df)
#       dG.pH8 <- soluble.orientation$best.soluble.score
#       pH8.df <- data.frame( pdb.id = pdb.id, pH.value = pH.value, dG = dG.pH8 )
#       main.df <- rbind( main.df, pH8.df )
#     }
#   }
# 
#   # Calculate the ddG of insertion values
#   data.df <- data.frame()
#   seq.tags <- unique(main.df$pdb.id)
#   for (seq in seq.tags) {
#     dG.pH4 <- main.df$dG[ which( main.df$pH.value == 4 & main.df$pdb.id == seq )]
#     dG.pH8 <- main.df$dG[ which( main.df$pH.value == 8 & main.df$pdb.id == seq )]
#     ddG.of.pH.shift = round( dG.pH4 - dG.pH8, 3)
#     variant.num <- unlist(strsplit(as.character(seq), "-|\\s")[1])[2]
#     exp.ddG <- exp.data$ddG[ which( exp.data$Sequence == seq ) ]
#     sub.df <- data.frame( seq = seq, var.no = variant.num, pred.ddG = ddG.of.pH.shift, exp.ddG = -exp.ddG )
#     data.df <- rbind( data.df, sub.df )
#   }
#   return(data.df)
# }
# 
# all.df <- read.mappings(workdir)
# write.table( all.df, file = paste( workdir, "processed_ddG_pH.txt", sep = "/") )

pH.df <- read.table( paste( workdir, "processed_ddG_pH.txt", sep = "/"), header = T )

sub.pH.df <- pH.df[ which( pH.df$var.no != "v1" & pH.df$var.no != "v14" & pH.df$var.no != "v16"), ]
print(lm(sub.pH.df$pred.ddG ~ sub.pH.df$exp.ddG ))
print(cor(-sub.pH.df$exp.ddG,sub.pH.df$pred.ddG))

p <- ggplot( data = pH.df, aes( x = -exp.ddG, y = pred.ddG, label = var.no ) ) +
  theme_bw() + 
  background_grid() +
  geom_abline( slope = 1.338, intercept = -11.15, color = "gray50", size = 0.35 ) + 
  geom_point( size = 0.75 ) +
  geom_text_repel( size = 2.5, force = 3 ) + 
  scale_x_continuous( "Experiment (kcal/mol)", limits = c(0, 3.0), expand = c(0,0) ) + 
  scale_y_continuous( "Predicted (REU)", limits = c(-20, -5), expand = c(0,0) ) + 
  theme( text = element_text( size = 7, color = "black" ), 
         axis.text = element_text( size = 7, color = "black" ), 
         axis.line = element_blank(), 
         axis.ticks = element_line( size = 0.35, color = "black" ), 
         panel.border = element_rect( color = "black" ))
print(p)
save_plot( paste( workdir, "C5_pH_helical_peptide.pdf", sep = "/"), p, units = "in", base_width = 2.24, base_height = 1.8 )

