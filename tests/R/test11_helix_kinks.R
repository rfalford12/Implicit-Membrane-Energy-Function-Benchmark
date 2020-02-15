#' Analyze helix kink angle prediction (Test #11)
#' 
#' Computes analysis for helix-kink angle prediction including: 
#'    - Kink angle vs. score plots
#'
#' @author Rebecca Alford <ralford3[at]jhu.edu>
#' 

library(cowplot)
library(reshape2)
library(ggplot2)

#' Read kink angle and score data into a dataframe
#'
#' @param score.file Path to scorefile including energies for each model
#' @param kinks.file Path to CSV file including the kink angle per helix and per model
#' @param is_gpcr Boolean variable that indicates if the target is a GPCR
#' @param is_exception Boolean variable that indicates the target uses a non-A chain
#'
#' @return Dataframe with scores and kink angles for each model of the target
#'
#' @examples
#' read.kink.angle.and.score.data( "3f7v_f19.sc", "3f7v_angles.csv", FALSE, FALSE )
#' 
read.kink.angle.and.score.data <- function( score.file, kinks.file, is_gpcr, is_exception ) { 
  
  scores.df <- read.table( paste( workdir, score.file, sep = "/"), header = T )
  kinks.df <- read.table( paste( workdir, kinks.file, sep = "/"), header = T )
  
  # Create a master data frame
  df <- data.frame()
  
  # Get score data
  scores <- scores.df$total_score
  
  # Will make dataframes for the most relevant angles (either 3/4 or 5/6/7)
  if ( !is_gpcr ) {
    
    helix.1.df <- kinks.df[ which( kinks.df$helix_num == 1 ), ]
    helix.2.df <- kinks.df[ which( kinks.df$helix_num == 2 ), ]
    if ( is_exception ) { 
      helix.2.df <- kinks.df[ which( kinks.df$helix_num == 3 ), ]
    }
    
    # sort by decoy number, ascending
    helix.1.df.sorted <- helix.1.df[ order( helix.1.df$model_number ), ]
    helix.2.df.sorted <- helix.2.df[ order( helix.2.df$model_number ), ]
    
    # assign scores column
    helix.1.df.sorted$total_score <- scores
    helix.2.df.sorted$total_score <- scores
    
    df <- rbind( helix.1.df.sorted, helix.2.df.sorted )
    
  } else {
    
    helix.5.df <- kinks.df[ which( kinks.df$helix_num == 5 ), ]
    helix.6.df <- kinks.df[ which( kinks.df$helix_num == 6 ), ]
    helix.7.df <- kinks.df[ which( kinks.df$helix_num == 7 ), ]
    
    # sort by decoy number, ascending
    helix.5.df.sorted <- helix.5.df[ order( helix.5.df$model_number ), ]
    helix.6.df.sorted <- helix.6.df[ order( helix.6.df$model_number ), ]
    helix.7.df.sorted <- helix.7.df[ order( helix.7.df$model_number ), ]
    
    helix.5.df.sorted$total_score <- scores
    helix.6.df.sorted$total_score <- scores
    helix.7.df.sorted$total_score <- scores 
    
    df <- rbind( helix.5.df.sorted, helix.6.df.sorted, helix.7.df.sorted )
    
  }
  return(df)
}

#' Make a plot of kink angle vs. score for a given target
#'
#' @param is_gpcr Boolean variable that indicates if the target is a GPCR
#' @param is_exception Boolean variable that indicates the target uses a non-A chain
#' @param native.angle.1 Native kink angle for first helix (1 or 5)
#' @param native.angle.2 Native kink angle for second helix (2 or 6)
#' @param native.angle.3 Native kink angle for third helix (3 or 7)
#' @param angle.limits Range for X-axis showing angles
#'
#' @return Plot object
#'
#' @examples
#' make.helix.kink.plot(df, TRUE, FALSE, 11.605, 35.813, 33.205, c(0,90) )
#' 
make.helix.kink.plot <- function( df, is_gpcr, is_exception, native.angle.1, native.angle.2, native.angle.3, angle.limits) {
  
  if ( !is_gpcr ) {
    if ( is_exception ) {
      p <- ggplot() + 
        theme_bw() + 
        background_grid() + 
        geom_point( data = df[ which( df$helix_num == 1 ), ], aes( x = kink_angle, y = total_score ), color = "#e41a1c", size = 0.25 ) + 
        geom_point( data = df[ which( df$helix_num == 3 ), ], aes( x = kink_angle, y = total_score ), color = "#253494", size = 0.25 ) + 
        geom_vline( xintercept = native.angle.1, color = "#e41a1c", linetype = "dashed" ) + 
        geom_vline( xintercept = native.angle.2, color = "#253494", linetype = "dashed" ) + 
        scale_x_continuous( "Kink Angle (degrees)",expand = c(0,0), limits = angle.limits ) + 
        scale_y_continuous("Score (REU)" ) + 
        theme( text = element_text( size = 8 ), 
               axis.text = element_text( size = 8 ))
      return(p)
    }
    p <- ggplot() + 
      theme_bw() + 
      background_grid() + 
      geom_point( data = df[ which( df$helix_num == 1 ), ], aes( x = kink_angle, y = total_score ), color = "#e41a1c", size = 0.25 ) + 
      geom_point( data = df[ which( df$helix_num == 2 ), ], aes( x = kink_angle, y = total_score ), color = "#253494", size = 0.25 ) + 
      geom_vline( xintercept = native.angle.1, color = "#e41a1c", linetype = "dashed" ) + 
      geom_vline( xintercept = native.angle.2, color = "#253494", linetype = "dashed" ) + 
      scale_x_continuous( "Kink Angle (degrees)",expand = c(0,0), limits = angle.limits ) + 
      scale_y_continuous("Score (REU)" ) + 
      theme( text = element_text( size = 8 ), 
             axis.text = element_text( size = 8 ))
    return(p)
  }
  
  p <- ggplot() + 
    theme_bw() + 
    background_grid() + 
    geom_point( data = df[ which( df$helix_num == 5 ), ], aes( x = kink_angle, y = total_score ), color = "#e41a1c", size = 0.25 ) + 
    geom_point( data = df[ which( df$helix_num == 6 ), ], aes( x = kink_angle, y = total_score ), color = "#253494", size = 0.25 ) + 
    geom_point( data = df[ which( df$helix_num == 7 ), ], aes( x = kink_angle, y = total_score ), color = "#984ea3", size = 0.25 ) + 
    geom_vline( xintercept = native.angle.1, color = "#e41a1c", linetype = "dashed" ) + 
    geom_vline( xintercept = native.angle.2, color = "#253494", linetype = "dashed" ) + 
    geom_vline( xintercept = native.angle.3, color = "#984ea3", linetype = "dashed" ) + 
    scale_x_continuous( "Kink Angle (degrees)",expand = c(0,0), limits = angle.limits ) + 
    scale_y_continuous("Score (REU)" ) + 
    theme( text = element_text( size = 8 ), 
           axis.text = element_text( size = 8 ))
  return(p)
  
  
}
