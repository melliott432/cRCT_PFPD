
library(R2jags)
library(rjags)
library(gridExtra)
library(truncnorm)
library(tidyverse)
library(grid)
library(ggplot2)
library(beepr)

setwd("D:/UCSD/Thesis/Bayesian cRCT design")

source("code/function_Estimate.R", echo=TRUE)
source("code/function_GetPowerDistn.R", echo=TRUE)
source("code/function_CalcPower.R", echo=TRUE)

################################################################################
### Format Myanmar data

# Trial data frame columns: 
# "arm_id" ----------------- binary  # 1 = control, 2 = treatment
# "village_id" ------------- integer # 1, 2, 3, 4, 5, 6, ...
# "adherence_final_indic" -- binary  # 0 = non-adherence, 1 = adherence
# "followup_final_indic" --- binary  # 0 = missing, 1 = observed

load("data/data_Bangladesh_pwr.RData")

df <- Bangladesh_pwr.df
df$village_id <- as.numeric(df$village_id)
df <- df %>%
  arrange(village_id)

################################################################################
### Run estimation

run_name = "myanmar"
estimates <- Estimate(run_name = run_name, df = df)
estimates 

b0 <- mean(estimates[,"beta0"])
b1 <- mean(estimates[,"beta1"])

save(estimates, file = "results/application/Bangladesh_estimates.RData")
################################################################################
### Power calculations

# Set operating characteristics and param estimates

alpha = 0.05   # type I error rate
n_ind = 15     # number of individuals per cluster
n_village = 16  # number of clusters per arm

parms = list("alpha" = alpha,                              
             "n_ind" = n_ind ,                             
             "n_village" = n_village)         

b0 <- median(estimates[,"beta0"])
b1 <- median(estimates[,"beta1"])

# Calculate Bayesian push-forward power distribution
power <- GetPower(run_name = "Bangladesh", estimates = estimates, parms = parms, B = 5000)

power.bayes <- power$power

power.freq <- CalcPower(alpha = 0.05, n = n_village, m = n_ind, 
                        pii = median(estimates[,"pii"]), 
                        tau = median(estimates[,"tau"]), 
                        delta = median(estimates[,"delta"]),
                        rho = median(estimates[,"rho_y"]), 
                        p1 = exp(b0)/(1+exp(b0)), 
                        p2 = exp(b0+b1)/(1+exp(b0+b1)))

################################################################################

# Save results

save(power.bayes, file = "results/application/Bangladesh_power_bayes.RData")
save(power.freq, file = "results/application/Bangladesh_power_freq.RData")
