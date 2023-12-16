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
# Main
###

# Read input
activity <- vroom(ACTIVITY, show_col_types = FALSE, col_names=TRUE)

# Number of CHIP files
nchip <- nrow(vroom(NCHIP, delim="\t", col_names=TRUE, show_col_types = FALSE))

# Main plot
 p1 <- ggplot(activity, aes(x=percentage, y=activity_winsor)) +
    geom_smooth(colour="blue", se=FALSE)

# Left CI
p2 <- ggplot(activity, aes(x=percentage, y=ci95_lbound_winsor)) + 
    geom_smooth(colour="black", se=FALSE) 

# Right CI
p3 <- ggplot(activity, aes(x=percentage, y=ci95_rbound_winsor)) + 
    geom_smooth(colour="black", se=FALSE)

# Extract smoothed lines from CI plots
main <- ggplot_build(p1)[[1]]
ci_l <- ggplot_build(p2)[[1]]
ci_r <- ggplot_build(p3)[[1]]

# Make plot data
plot_data <-data.frame(x = main[[1]]$x,
                  ymin = ci_l[[1]]$y,
                  ymax = ci_r[[1]]$y, 
                  mainy = main[[1]]$y)

# Final
b <- ggplot(plot_data, aes(x=x,  y=mainy, ymin = ymin, ymax = ymax)) +
    geom_line(colour="#4DBBD5FF", linewidth=1) +
    geom_point(alpha=0.2) + 
    geom_ribbon(fill="grey", alpha=0.4) +
    xlab("PWM Score Percentile") + 
    ylab("Proportion Unibind Support") +
    xlim(79.99, 90.01) +
    ylim(0, 1) +
    theme_minimal(base_size = 20)

# Pub params
b <- ggpar(b, palette="npg", title=PROFILE, subtitle=nchip) + 
    theme_pubr(base_size=16, base_family="sans")

# Export
ggarrange(b) %>%
  ggexport(filename = OUTPUT)