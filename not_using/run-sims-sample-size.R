# Simulations for push forward power with varying sample sizes

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

source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_SimData.R")
source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_Estimate.R")
source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_GetPower.R")

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code/sims_sample_size")

# Parameters that don't vary
alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
pii = 0.75 # follow-up rate
tau = 0.25 # ICC for missingness
delta = 0.2 # treatment effect
sd_mult = 0.01 # standard deviation for outcome model

# Sample size combinations to try
n_clusters = c(6, 10, 16, 20)
cluster_size = c(25, 50, 100, 200)
#n_clusters <- ifelse(n_clusters %% 2 == 0, n_clusters, n_clusters+1)
sim_sample_size <- expand.grid(n_clusters = n_clusters, cluster_size = cluster_size)
#sim_sample_size$total_samp_size <- sim_sample_size$n_clusters * sim_sample_size$cluster_size
#sim_sample_size <- sim_sample_size |> filter(total_samp_size < 5000) 

# Data frame to save info about each simulation
labels <- c("run_name", "n_village", "n_ind", "n_total", "sim_follow_up_rate", "sim_attrition_icc", "sim_treatment_eff", "run.time")
sim_info <- as.data.frame(matrix(data = NA, nrow = nrow(sim_sample_size), ncol = length(labels)))
colnames(sim_info) <- labels

#### Only change input above this line ######################################### 
################################################################################
################################################################################

# Iterate through combos
for(i in 1:nrow(sim_sample_size)){
  
  # Check if combo has been run
  if(is.na(sim_info$run_name[i])){
    
    start.time <- Sys.time() # track run time
    
    ############################################################################
    #### SIMULATE DATA #########################################################
    
    n_village <- sim_sample_size$n_clusters[i]
    n_ind <- sim_sample_size$cluster_size[i]
    n_ind_total <- n_village*n_ind
    
    run_name <- paste(Sys.Date(), "_", n_village,"_", n_ind, sep = "") # label for this run
    
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
    
    sim_info[i,] <- c(run_name, n_village, n_ind, n_ind_total, round(c(pi.sim, tau.sim, delta.sim, run.time),3))

  }
}

save(sim_info, file = paste(Sys.Date(),"info.RData",sep = ""))

#################################################################################

beep(1)





