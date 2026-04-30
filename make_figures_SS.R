##### Create figures of SAMPLE SIZE results ####################################

library(grobblR)
library(gridExtra)
library(dplyr)
library(png)
library(viridis)
library(tidyverse)
library(ggplot2)


setwd("D:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

load("sample_size_info.RData")
load("sample_size_power_df.RData")

################################################################################
### Figure X: Power Distributions ##############################################
################################################################################

colors <- viridis(4)

power.df$n_village <- as.factor(power.df$n_village)
power.df$n_ind <- as.factor(power.df$n_ind)

png("sample_size_figure.png")
ggplot(power.df, aes(x = avg_median_power, colour = n_ind)) +
  facet_wrap(~ n_village) +
    geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75) +
    labs(title = paste("Kernel Density Plot of Power, 6 clusters", sep = ""),
         x = "Power",
         y = "Density") +
    xlim(0,1) + 
    theme_light() +
    scale_color_manual(values = colors)
dev.off()
