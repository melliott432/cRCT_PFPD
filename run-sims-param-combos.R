# Simulations for push forward power with varying parameter values

################################################################################
##### SET UP ###################################################################

library(R2jags)
library(rjags)
library(gridExtra)
library(truncnorm)
library(tidyverse)
library(grid)
library(ggplotify)
library(beepr)

source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_SimData.R", echo=TRUE)
source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_Estimate.R", echo=TRUE)
source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_GetPower.R", echo=TRUE)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code/sims_param_combos")

### CONSTANT OPERTAING CHARACTERISTICS #########################################
alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
sd_mult = 0.05 # standard deviation for outcome model
n_village = 10 
n_ind = 200

### PARAMS TO VARY #############################################################
pii = c(0.50, 0.75, 0.90)
tau = c(0.25, 0.50, 0.75)
delta = c(0.1, 0.2, 0.3)
sim_params <- expand.grid(pii = pii, tau = tau, delta = delta)

### SAVE RUN INFO ##############################################################
labels <- c("run_name", "follow_up_rate", "attrition_icc", "treatment_eff", "sim_follow_up_rate", "sim_attrition_icc", "sim_treatment_eff", "run.time")
sim_info <- as.data.frame(matrix(data = NA, nrow = nrow(sim_sample_size), ncol = length(labels)))
colnames(sim_info) <- labels

################################################################################

# Iterate through combos
for(i in 1:nrow(sim_params)){
  
  start.time <- Sys.time() # track run time
    
  ### Set Parameters ###
  
  pii <- sim_params[i, "pii"] # follow-up rate
  tau <- sim_params[i, "tau"] # attrition ICC
  delta <- sim_params[i, "delta"] # treatment effect
    
  parms <- list("n_village" = n_village, 
                  "n_ind" = n_ind, 
                  "n_arm" = n_arm, 
                  "pii" = pii, 
                  "tau" = tau, 
                  "delta" = delta, 
                  "sd_mult" = sd_mult)
    

  ### Run simulation and estimation ###
  
  run_name <- paste(Sys.Date(), pii, tau, delta, sep = "_") 

  simdata <- SimData(run_name = run_name, parms = parms)
    
  estimates <- SimEstimate(run_name = run_name, simdata = simdata)
    
  power <- GetPower(run_name = run_name, simdata = simdata, estimates = estimates)
  
  png(filename = paste("plots/power/", run_name, "_power", ".png", sep = ""))
  plot(density(power), main = paste("n = ", n_village, ", m = ", n_ind, sep = ""), xlab = "Power")
  dev.off()
  
  ### Check Sim Data ### 
  
  pi.sim <- mean(simdata$followup_final_indic)
    
  pr = simdata |>
    group_by(village_id) |>
    summarize(prop = sum(followup_final_indic)/length(followup_final_indic))
    
  tau.sim <- var(pr$prop)/(pi.sim*(1-pi.sim))
    
  a.sim <- (pi.sim*(1-tau.sim))/tau.sim
  b.sim <- (1/tau.sim) - 1 - ((pi.sim*(1-tau.sim))/tau.sim)
    
  p1 <- mean(simdata$adherence_final_indic[which(simdata$arm_id==1)])
  p2 <- mean(simdata$adherence_final_indic[which(simdata$arm_id==2)])
  delta.sim <- p2 - p1
  

  # Get run time
  end.time <- Sys.time()
  run.time <- end.time - start.time
    
  # Save info about this run
  sim_info[i,] <- c(run_name, pii, tau, delta, round(c(pi.sim, tau.sim, delta.sim, run.time),3))

}

save(sim_info, file = "info.RData")

#################################################################################







