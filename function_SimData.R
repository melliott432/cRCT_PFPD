### Simulates CRT data (adherence and attrition) for one trial

SimData <- function(run_name, parms){
  
  require(geeCRT)
  require(dplyr)
  require(logitnorm)
  require(mvtnorm)
  
  n_village = parms$n_village
  n_ind = parms$n_ind
  n_arm = parms$n_arm
  pii = parms$pii
  tau = parms$tau
  delta = parms$delta
  sd_mult = parms$sd_mult

  n_ind_total = n_village*n_ind
  
  beta = logit(delta + 0.5)
  
  seed = as.numeric(pii*100 + tau*100 + delta*100) + n_village + n_ind
  set.seed(seed)
    
  ##############################################################################
  ## ATTRITION
  
  # Beta distribution for cluster-level proportions
  a = (pii*(1-tau))/tau
  b = (1/tau) - 1 - a
  pr_cluster = rbeta(n = n_village, shape1 = a, shape2 = b)
  village_flwup_rate = rep(pr_cluster, each = n_ind)
    
  # Binomial distribution for individual missingness
  followup_final_indic <- c()
  for(i in 1:n_village){
    Rvec <- rbinom(n = n_ind, size = 1, prob = pr_cluster[i])
    followup_final_indic <- c(followup_final_indic, Rvec)
  }
  ## R = 1 is observed, R = 0 is missing

  ##############################################################################
  ## ADHERENCE

  # Create std devs for random effects
  village_sd = 1*sd_mult # village-level RE
  clinic_sd = 1*sd_mult # clinic-level RE
  household_sd = 1*sd_mult # household-level RE
  individual_sd = 1*sd_mult # individual-level RE
  
  # Create ID vectors
  ind_id = c(1:n_ind_total)
  arm_id = c(rep(c(1:n_arm), each = n_ind_total/n_arm))
  village_id = c(rep(c(1:n_village), each = n_ind))
  
  # Intervention (arm) effect
  arm_eff = c(rep(0, length = n_ind_total/2), rep(beta, length = n_ind_total/2))
  
  # Village (cluster) effect
  village_eff = rnorm(n_village, mean = 0, sd = village_sd)
  village_eff = rep(village_eff, each = n_ind)
  
  # Individual effect
  individual_eff = rnorm(n_ind_total, mean = 0, sd = individual_sd)
  
  # Sample baseline adherence indicators
  adherence_base_logit = village_eff + individual_eff
  adherence_base_prob = invlogit(adherence_base_logit)
  adherence_base_indic = rbinom(n_ind_total, 1, prob = adherence_base_prob)
  
  # Sample follow-up adherence indicators
  adherence_final_logit = village_eff + individual_eff + arm_eff
  adherence_final_prob = invlogit(adherence_final_logit)
  adherence_final_indic = rbinom(n_ind_total, 1, prob = adherence_final_prob)
  
  adherence_final_indic[which(followup_final_indic == 1)] = NA
  
  ##############################################################################
  # COMBINE AND RETURN
  
  sim_trial_df <- data.frame(ind_id, village_id, arm_id, 
                             individual_eff, village_eff, arm_eff, 
                             adherence_base_logit, adherence_base_prob, adherence_base_indic, 
                             adherence_final_logit, adherence_final_prob, adherence_final_indic,
                             village_flwup_rate, followup_final_indic)
  
  filename <- paste("datasets/simdata_", run_name, ".RData", sep="")
  save(sim_trial_df, file = filename)
  
  return(sim_trial_df)

}

