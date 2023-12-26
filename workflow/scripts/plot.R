library(repr)
library(dplyr)
library(vroom)
library(ggplot2)
library(pracma)
library(patchwork)
library(ggpubr)
library(ggsci)
library(ggExtra)


# Snakemake
ACTIVITY <- snakemake@input[["activity"]]
PROFILE <- snakemake@params[["profile"]]
NCHIP <- snakemake@input[["nchip"]]
OUTPUT <- snakemake@output[[1]]

###
# Function
###

plot_activity <- function(dataset, xvar, mean_var, lbound_var, rbound_var, plot_xlab, plot_ylab, plot_title, num_chip) {
    # Main
    main <- ggplot(dataset, aes_string(x=xvar, y=mean_var)) +
        geom_smooth(colour="blue", se=FALSE)

    # Left CI
    lbnd <- ggplot(dataset, aes_string(x=xvar, y=lbound_var)) + 
        geom_smooth(colour="black", se=FALSE) 

    # Right CI
    rbnd <- ggplot(dataset, aes_string(x=xvar, y=rbound_var)) + 
        geom_smooth(colour="black", se=FALSE)

    # Extract smoothed lines from CI plots
    main_data <- ggplot_build(main)[[1]]
    lbnd_data <- ggplot_build(lbnd)[[1]]
    rbnd_data <- ggplot_build(rbnd)[[1]]

    # Construct into matrix
    plot_data <-data.frame(x = main_data[[1]]$x,
                  ymin = lbnd_data[[1]]$y,
                  ymax = rbnd_data[[1]]$y, 
                  mainy = main_data[[1]]$y)

    # Plot together
    plot <- ggplot(plot_data, aes(x=x,  y=mainy, ymin = ymin, ymax = ymax)) +
        geom_line(colour="#4DBBD5FF", linewidth=1) +
        geom_point(alpha=0.2) + 
        geom_ribbon(fill="grey", alpha=0.4) +
        xlab(plot_xlab) + 
        ylab(plot_ylab) +
        xlim(79.99, 100.01) +
        ylim(0, 1) +
        theme_minimal(base_size = 20)

    # Pub params
    plot <- ggpar(plot, palette="npg", title=plot_title, subtitle=num_chip) + 
        theme_pubr(base_size=16, base_family="sans")
        
    # Return final plot
    return (plot)
}

###
# Main
###

# Read input
activity <- vroom(ACTIVITY, show_col_types = FALSE, col_names=TRUE)

# Number of CHIP files
nchip <- nrow(vroom(NCHIP, delim="\t", col_names=TRUE, show_col_types = FALSE))

# Make plots
a <- plot_activity(dataset = activity, xvar = "percentage", mean_var = "activity",
 lbound_var = "ci95_lbound", rbound_var = "ci95_rbound", plot_xlab = "PWM Score Percentile", 
 plot_ylab = "Proportion Unibind Support", plot_title = PROFILE, num_chip = nchip)


# b <- plot_activity(dataset = activity, xvar = "percentage", mean_var = "activity_winsor",
#  lbound_var = "ci95_lbound", rbound_var = "ci95_rbound", plot_xlab = "PWM Score Percentile", 
#  plot_ylab = "Proportion Unibind Support", plot_title = "WINSOR", num_chip = nchip)

# Combine plots
options(repr.plot.width=24, repr.plot.height=16)
c <- a #+ b

# Export
ggarrange(c) %>%
  ggexport(filename = OUTPUT)