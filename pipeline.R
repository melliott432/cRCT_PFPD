################################################################################
#### Environment set-up

setwd("D:/UCSD/Thesis/Bayesian cRCT design")

library(gridExtra)
library(tidyverse)

source("code/function_SimulateData.R", echo=TRUE)
source("code/function_Estimate.R", echo=TRUE)
source("code/function_GetPowerDistn.R", echo=TRUE)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

################################################################################
#### Simulation settings to run

# Load in info of simulations to run
if(TYPE == "PC"){
  load("param_combos_info.RData")
}else if(TYPE == "SS"){
  load("sample_size_info.RData")
}else{
  print("Specify which simulations to run (PC or SS)")
}

siminfo <- mutate(siminfo, estimate_file_name = paste("posterior_", run_name, ".RData", sep = ""),
                  simdata_file_name = paste("simdata_", run_name, ".RData", sep = ""))

files.done <- list.files()

torun <- which(siminfo$estimate_file_name %in% files.done == FALSE)

################################################################################
#### Iterate through parameter combinations

for(j in 1:length(torun)){
  
  #start.time <- Sys.time() # track run time

  run_name <- siminfo$run_name[torun[j]] # set name for this run
    
  parms <- list("n_village" = siminfo$n_village[torun[j]], 
                "n_ind" = siminfo$n_ind[torun[j]], 
                "n_arm" = siminfo$n_arm[torun[j]], 
                "sd_mult" = siminfo$sd_mult[torun[j]],
                "alpha" = siminfo$alpha[torun[j]],
                "pii" = siminfo$pii[torun[j]], 
                "tau" = siminfo$tau[torun[j]], 
                "delta" = siminfo$delta[torun[j]])
    
  #simdata <- SimData(run_name = run_name, parms = parms) # simulate data
  load(paste("simdata_", run_name, ".RData", sep = ""))
  simdata <- sim_trial_df
  
  #estimates <- Estimate(run_name = run_name, df = simdata) # estimate parameter posteriors
  load(paste("posterior_", run_name, ".RData", sep = ""))
  estimates <- posterior
  
  power <- GetPower(run_name = run_name, df = simdata, estimates = estimates, parms = parms) # calculate power push-forward

  #siminfo[torun[j],"run_time"] <- round(Sys.time() - start.time,3) # calculate total run time for this rep
    
  #save(siminfo, file = paste("param_combos_", "info.RData", sep = "")) # update sim info and save
}
 