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
METADATA <- snakemake@input[[1]]
 <- snakemake@input[[2]]
OUTPUT <- snakemake@output[[1]]

###
# Main
###

# Read input
metadata <- vroom(METADATA, show_col_types = FALSE, col_names=TRUE)