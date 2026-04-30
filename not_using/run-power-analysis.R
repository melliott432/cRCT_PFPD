
library(R2jags)
library(rjags)
library(gridExtra)
library(truncnorm)
library(tidyverse)
library(grid)
library(ggplot2)
library(beepr)

source("function_SimData.R", echo=TRUE)
source("function_Estimate.R", echo=TRUE)
source("function_GetPower.R", echo=TRUE)

run_name = "test_run"

################################################################################

## Arguments
    # df - data frame of trial data
    # needs the following columns, all numeric:
      # "arm_id" ----------------- binary
      # "village_id" ------------- integer
      # "adherence_base_indic" -- binary
      # "adherence_final_indic" -- binary
      # "followup_final_indic" --- binary 

estimates <- SimEstimate(run_name = run_name, df = df)


################################################################################

## Arguments   
    # parms - named list:
      # alpha = parms$alpha
      # n_ind = parms$n_ind
      # n_village = parms$n_village
      # pii = parms$pii
      # tau = parms$tau
      # delta = parms$delta
      # pii, tau, delta can be point estimates from above

parms = list(alpha = 0.05,
             n_ind = nrow(df),
             n_village = length(unique(df$village_id)),
             pii = mean(estimates[,"pii"]),
             tau = mean(estimates[,"tau"]),
             delta = mean(estimates[,"delta"]))

power <- GetPower(run_name = run_name, df = df, estimates = estimates, parms = parms)
        
################################################################################

plot_name = "myplot"
  
power.df <- data.frame("power" = power[["power.dist"]]) 

ggplot(power.df, aes(x = power)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75) +
  labs(title = "Kernel Density Plot of Power",
       x = "Power", y = "Density") +
  xlim(0,1) + 
  theme_light()

ggsave(filename = paste(plot_name, ".png", sep =""))










