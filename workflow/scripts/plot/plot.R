# Load required libraries
library(repr)
library(dplyr)
library(vroom)
library(ggplot2)
library(pracma)
library(patchwork)
library(ggpubr)
library(ggsci)
library(ggExtra)
library(DescTools)

# Snakemake
ACTIVITY <- snakemake@input[["activity"]]
PROFILE <- snakemake@params[["profile"]]
DATASET <- snakemake@input[["dataset"]]
PLOT_OUTPUT <- snakemake@output[["plt"]]
AUC_OUTPUT <- snakemake@output[["auc"]]

# ------------- #
# Functions     #
# ------------- #

plot_activity <- function(dataset, xvar, mean_var, lbound_var, rbound_var, plot_xlab, plot_ylab, plot_title, dataset_id) {
    tryCatch(
         expr = {
            options(repr.plot.width=24, repr.plot.height=16)
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

            # Grab auc
            auc <- AUC(main_data[[1]]$x, main_data[[1]]$y)

            # Save auc
            write(as.character(auc), AUC_OUTPUT)

            # Construct into matrix
            plot_data <-data.frame(x = main_data[[1]]$x,
                        ymin = lbnd_data[[1]]$y,
                        ymax = rbnd_data[[1]]$y, 
                        mainy = main_data[[1]]$y)

            min(plot_data$x)
            # Plot together
            plot <- ggplot(plot_data, aes(x=x,  y=mainy, ymin = ymin, ymax = ymax)) +
                geom_line(colour="#4DBBD5FF", linewidth=1) +
                geom_point(alpha=0.2) + 
                geom_ribbon(fill="grey", alpha=0.4) +
                xlab(plot_xlab) + 
                ylab(plot_ylab) +
                xlim(min(plot_data$x)-0.01, 100.01) +
                ylim(0, 1) +
                theme_minimal(base_size = 20)

            # Pub params
            plot <- ggpar(plot, palette="npg", title=plot_title, subtitle=dataset_id) + 
                theme_pubr(base_size=16, base_family="sans")
                
        },
        error = function(e){
            message('Caught an error!')
            print(e)
            plot <- ggplot(activity, aes(percentage, activity)) + geom_blank()
        },
       finally = {
        message('Saving plot...')
            ggarrange(plot) %>%
                 ggexport(filename = PLOT_OUTPUT)
        }
    )    
}



# ------------- #
# Main          #
# ------------- #

# Read input
activity <- vroom(ACTIVITY, show_col_types = FALSE, col_names=TRUE)

# Make plots
plot_activity(dataset = activity, xvar = "perc", mean_var = "activity",
                lbound_var = "ci95_lbound", rbound_var = "ci95_rbound", plot_xlab = "Percent Max PWM Score", 
                plot_ylab = "Proportion Unibind Support", plot_title = PROFILE, dataset_id = DATASET)
