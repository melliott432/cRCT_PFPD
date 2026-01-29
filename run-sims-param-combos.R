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

setwd("E:/UCSD/Thesis/Bayesian cRCT design/code")

source("function_SimData.R", echo=TRUE)
source("function_Estimate.R", echo=TRUE)
source("function_GetPower.R", echo=TRUE)

setwd("sims_param_combos")

### CONSTANT OPERTAING CHARACTERISTICS #########################################
alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
sd_mult = 0.05 # standard deviation for outcome model
n_village = 10 
n_ind = 50

### PARAMS TO VARY #############################################################
pii = seq(0.1, 0.95, by = 0.05)
tau = seq(0.1, 0.95, by = 0.05)
delta = seq(0.1, 0.95, by = 0.05)
sim_params <- expand.grid(pii = pii, tau = tau, delta = delta)

### RUN INFO ##############################################################
labels <- c("run_name", "follow_up_rate", "attrition_icc", "treatment_eff", "sim_follow_up_rate", "sim_attrition_icc", "sim_treatment_eff", "run.time")
sim_info <- as.data.frame(matrix(data = NA, nrow = nrow(sim_params), ncol = length(labels)))
colnames(sim_info) <- labels

done = list.files(path = "power")

sim_info = sim_params |> 
  dplyr::mutate(
    run_name <- paste("2026-01-28", pii*100, tau*100, delta*100, "n", n_village, "m", n_ind, sep = "_")
  ) |> 
  filter(
  #  !(FILE %in% done)
  #) 

#### Only change input above this line ######################################### 
################################################################################
################################################################################

# Iterate through combos
for(i in 1:nrow(sim_params)){
  
  # Check if combo has been run
  if(is.na(sim_info$run_name[i])){
    
    start.time <- Sys.time() # track run time
    
    ############################################################################
    #### SIMULATE DATA #########################################################
    
    pii <- sim_params$pii[i]
    tau <- sim_params$tau[i]
    delta <- sim_params$delta[i]
    
    parms <- list("n_village" = n_village, 
                  "n_ind" = n_ind, 
                  "n_arm" = n_arm, 
                  "pii" = pii, 
                  "tau" = tau, 
                  "delta" = delta, 
                  "sd_mult" = sd_mult,
                  "alpha" = alpha)
    
    simdata <- SimData(run_name = run_name, parms = parms)
    
    estimates <- SimEstimate(run_name = run_name, simdata = simdata)
    
    power <- GetPower(run_name = run_name, simdata = simdata, estimates = estimates, parms = parms)
    
    
    # Check pi
    pi.sim <- mean(simdata$followup_final_indic)
    
    # Check tau
    pr = simdata |>
      group_by(village_id) |>
      summarize(prop = sum(followup_final_indic)/length(followup_final_indic))
    
    tau.sim <- var(pr$prop)/(pi.sim*(1-pi.sim))
    
    a.sim <- (pi.sim*(1-tau.sim))/tau.sim
    b.sim <- (1/tau.sim) - 1 - ((pi.sim*(1-tau.sim))/tau.sim)
    
    # Check treatment effect
    p1 <- mean(simdata$adherence_final_indic[which(simdata$arm_id==1)])
    p2 <- mean(simdata$adherence_final_indic[which(simdata$arm_id==2)])
    delta.sim <- p2 - p1
    
    end.time <- Sys.time()
    
    run.time <- end.time - start.time
    
    sim_info[i,] <- c(run_name, pii, tau, delta, round(c(pi.sim, tau.sim, delta.sim, run.time),3))
    
  }
}

save(sim_info, file = paste(Sys.Date(), "_n", n_village, "_m", n_ind, "_", "info.RData", sep = ""))
 
