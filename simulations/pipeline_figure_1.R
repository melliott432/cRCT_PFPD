
library(gridExtra)
library(tidyverse)

setwd("E:/UCSD/Thesis/Bayesian cRCT design/code")

source("function_SimulateData.R", echo=TRUE)
source("function_Estimate.R", echo=TRUE)
source("function_GetPowerDistn.R", echo=TRUE)
source("function_CalcPower.R", echo = TRUE)


n_reps = 5 # number of replicates to run

alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
sd_mult = 0.05 # standard deviation for outcome model
p1 = 0.3

################################################################################
### SAMPLE SIZE COMBINATIONS ###################################################

n_village = c(10, 12,14,18) # number of clusters total
n_ind = c(25, 50, 100, 200) # cluster size

pii = 0.42 # follow-up rate
tau = c(0.01, 0.06, 0.1, 0.3) # attrition ICC
delta = 0.2 # effect size


parms_grid <- expand.grid("rep" = 1:n_reps, "alpha" = alpha, "n_arm" = n_arm, 
                          "sd_mult" = sd_mult, "p1" = p1, 
                          "n_village" = n_village, "n_ind" = n_ind, 
                          "pii" = pii, "tau" =  tau, "delta" = delta)
                          
setwd("E:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

for(i in torun){
  
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
  
  print(paste("Done running row ", i, sep = " "))
}

################################################################################

power.df <- data.frame(parms_grid, "med_power" = rep(NA,nrow(parms_grid)), "freq_power" = rep(NA,nrow(parms_grid)))

for(i in 1:nrow(parms_grid)){
  
  run_name = paste(100*parms_grid[i,"pii"],100*abs(parms_grid[i,"tau"]),100*abs(parms_grid[i,"delta"]), parms_grid$rep[i], sep = "_")
  
  load(file = paste("power_", run_name, ".RData", sep = ""))
  
  power.df$med_power[i] <- median(power.dist$power)
  
  power.df$freq_power[i] <- CalcPower(alpha = parms_grid$alpha[i],
                                      n = parms_grid$n_village[i],
                                      m = parms_grid$n_ind[i],
                                      n_arms = parms_grid$n_arm[i],
                                      pii = parms_grid$pii[i],
                                      tau = parms_grid$tau[i],
                                      rho = mean(power.dist$rho),
                                      p1 = parms_grid$p1[i],
                                      p2 = parms_grid$p1[i] + parms_grid$delta[i])
  
}

power.df.1 <- power.df |>
  group_by(pii, tau, delta) |>
  summarise(avg_med_power = mean(med_power),
            avg_freq_power = mean(freq_power))

save(power.df.1, file = "power_results_summary.RData")
