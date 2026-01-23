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

source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_SimData.R", echo=TRUE)
source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_Estimate.R", echo=TRUE)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code/sims_sample_size")

# Parameters that don't vary
alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
pii = 0.75 # follow-up rate
tau = 0.25 # ICC for missingness
delta = 0.2 # treatment effect
sd_mult = 0.01 # standard deviation for outcome model

# Sample size combinations to try
n_clusters = seq(5, 200, by = 20)
n_clusters <- ifelse(n_clusters %% 2 == 0, n_clusters, n_clusters+1)
sim_sample_size <- expand.grid(n_clusters = n_clusters, cluster_size = round(seq(5, 160, length.out = 10)))
sim_sample_size$total_samp_size <- sim_sample_size$n_clusters * sim_sample_size$cluster_size
sim_sample_size <- sim_sample_size |> filter(total_samp_size < 5000) 

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
    
    run_name <- paste(Sys.Date(), n_village, n_ind, sep = "_") # label for this run
    
    # Call data sim function
    simdata <- SimData(n_village = n_village, n_ind = n_ind, n_arm = n_arm, pii=pii, tau = tau, delta = delta, sd_mult = sd_mult)
    filename <- paste("datasets/",paste(run_name, "RData", sep = "."), sep="")
    save(simdata, file = filename)
    
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
  
    ############################################################################
    #### ESTIMATE #######################################################
    
    estimates <- Sim_Estimate(simdata = simdata)
    
    ############################################################################
    #### CALCULATE POWER #######################################################
    
    pr = simdata |>
      group_by(village_id) |>
      summarize(prop = sum(adherence_base_indic)/length(adherence_base_indic))
    var.p = var(pr$prop)
    exp.p = mean(pr$prop)
    rho_y = var.p / (exp.p*(1-exp.p))
    
    p1 = mean(simdata$adherence_base_indic)
    p2 = delta + p1
    
    # Sample parameters from posterior distribution
    B = 1000
    a.B <- sample(estimates[,"a"], B)
    b.B <- sample(estimates[,"b"], B)
    delta.B <- sample(estimates[,"delta"],B)
    tau.B <- 1/(a.B+b.B+1)
    pii.B <- a.B/(a.B+b.B)
    
    G = 1 + (n_ind*pii.B - 1)*rho_y + (1-pii.B)*(1 + (n_ind-1)*tau.B)*rho_y
    
    term = n_village*n_ind*delta.B^2 *(pii.B/G) *(p1*(1-p1) + (delta.B+p1)*(1-(delta.B+p1)))^-1
    
    power <- pnorm(sqrt(term) - qnorm(1-alpha/2))
    

    # Power
    png(filename = paste("plots/power/", run_name, "_power", ".png", sep = ""))
    plot(density(power), main = paste("n = ", n_village, ", m = ", n_ind, sep = ""), xlab = "Power")
    dev.off()
    
    end.time <- Sys.time()
    
    run.time <- end.time - start.time
    
    sim_info[i,] <- c(run_name, n_village, n_ind, n_ind_total, round(c(pi.sim, tau.sim, delta.sim, run.time),3))

  }
}

save(sim_info, file = paste(Sys.Date(),"info.RData",sep = ""))

#################################################################################

beep(1)





