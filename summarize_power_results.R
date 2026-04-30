
library(viridis)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_CalcPower.R")

################################################################################
### PARAM COMBOS ###############################################################
################################################################################

load("param_combos_info.RData")

################################################################################
## Median of Bayesian push-forward #############################################

# Calculate median power for each rep
power_median <- rep(NA, nrow(siminfo))
for(i in 1:nrow(siminfo)){
  filename = paste("power_", siminfo$run_name[i], ".RData", sep = "")
  load(filename)
  power_median[i] <- median(power.bayes)
}
siminfo <- cbind(siminfo, power_median)

# Calculate average median power for each set 
power_avg_median <- siminfo |> 
  select(-run_name, -run_time) |>
  group_by(pii, tau, delta) |>
  summarize(avg_median_power = mean(power_median))
                      
################################################################################
## Frequentist point comparison ################################################

rhoy = (0.1^2) / ((0.1^2) + ((pi^2)/3))

power_freq <- siminfo |>
  select(-rep, -run_name, -run_time) |>
  distinct() |>
  mutate(freq_power = CalcPower(alpha = alpha, n = n_village, m = n_ind, 
                                pii = pii, tau = tau, delta = delta, rho = rhoy,
                                P0 = 0.3, P1 = (0.3+delta))) 

################################################################################
## Save power summary ##########################################################

power.summary <- merge(power_freq, power_avg_median)

save(power.summary, file = "param_combos_power_summary.RData")

################################################################################
################################################################################

################################################################################
### SAMPLE SIZE ################################################################
################################################################################

load("sample_size_info.RData")
source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_CalcPower.R", echo = TRUE)

################################################################################
## Median of Bayesian push-forward #############################################


power_median <- rep(NA, nrow(siminfo))

power.df <- matrix(NA, nrow=nrow(siminfo), ncol=1000)

for(i in 1:nrow(siminfo)){
  
  filename = paste("power_", siminfo$run_name[i], ".RData", sep = "")
  load(filename)
  
  power_median[i] <- median(power.bayes)
  
  power.df[i,] <- power.bayes
  
}

# Calculate average median power for each set 
power_avg_median <- siminfo |> 
  select(-run_name, -run_time) |>
  group_by(n_village, n_ind) |>
  summarize(avg_median_power = mean(power_median), lower = (mean(power_median) - 1.96*sd(power_median)), upper = mean(power_median) + 1.96*sd(power_median))

################################################################################
## Frequentist point comparison ################################################

rhoy = (0.1^2) / ((0.1^2) + ((pi^2)/3))

power_freq <- siminfo |>
  select(-rep, -run_name, -run_time) |>
  distinct() |>
  mutate(freq_power = CalcPower(alpha = alpha, n = n_village, m = n_ind, 
                                pii = pii, tau = tau, delta = delta, rho = rhoy,
                                P0 = 0.3, P1 = (0.3+delta))) 


################################################################################
## Save power summary ##########################################################

power.summary <- merge(power_freq, power_avg_median)
save(power.summary, file = "sample_size_power_summary.RData")

power.df <- cbind(siminfo, power.df)
save(power.df, file = "sample_size_power_df.RData")
