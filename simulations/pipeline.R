
library(gridExtra)
library(tidyverse)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code")

source("function_SimulateData.R", echo=TRUE)
source("function_Estimate.R", echo=TRUE)
source("function_GetPowerDistn.R", echo=TRUE)
source("function_CalcPower.R", echo = TRUE)

parms_grid <- expand.grid("rep" = 1:n_reps, "alpha" = alpha, "n_arm" = n_arm, 
                          "sd_mult" = sd_mult, "p1" = p1, 
                          "n_village" = n_village, "n_ind" = n_ind, 
                          "pii" = pii, "tau" =  tau, "delta" = delta)
                          
setwd("D:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

for(i in 1:nrow(parms_grid)){
  
  parms <- list("n_village" =  parms_grid$n_village[i], # total number of clusters
                "n_ind" = parms_grid$n_ind[i],     # common cluster size
                "n_arm"= parms_grid$n_arm[i],       # number of arms
                "pii" = parms_grid$pii[i],     # follow-up rate
                "tau" = parms_grid$tau[i],     # attrition ICC
                "delta" = parms_grid$delta[i],   # treatment effect (raw change in proportion)
                "sd_mult" = parms_grid$sd_mult[i], # common standard deviation,
                "p1" = parms_grid$p1[i], # baseline or control outcome proportion,
                "alpha" = 0.05
  )
  
  run_name = paste(parms_grid[i,"n_village"],parms_grid[i,"n_ind"],100*parms_grid[i,"pii"],100*abs(parms_grid[i,"tau"]),100*abs(parms_grid[i,"delta"]),parms_grid[i,"rep"],sep = "_")
  
  sim.df <- SimData(run_name = run_name, parms = parms)
  save(sim.df, file = paste("simdata_",run_name,".RData",sep=""))
  
  estimates <- Estimate(df = sim.df, run_name = run_name)
  save(estimates, file = paste("posterior_",run_name,".RData",sep=""))
  
  power.dist <- GetPower(run_name = run_name, parms = parms, estimates = estimates)
  save(power.dist, file = paste("power_", run_name, ".RData", sep = ""))
  
  print(paste("Done running", run_name, sep = " "))
}
