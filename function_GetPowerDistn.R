GetPower <- function(run_name, df, estimates, parms){
  
  source("D:/UCSD/Thesis/Bayesian cRCT design/code/function_CalcPower.R")
  
  # Define parameters
  alpha = parms$alpha # type I error rate aka significance level
  n_ind = parms$n_ind # cluster size
  n_village = parms$n_village # number of clusters
  
  # Number of MCMC iterations
  B = 1000

  # Sample from posterior
  pii.B <- sample(estimates[,"pii"], B)
  tau.B <- sample(estimates[,"tau"], B)
  delta.B <- sample(estimates[,"delta"],B)
  rho_y.B <- sample(estimates[,"rho_y"],B)
  beta0.B <- sample(estimates[,"beta0"],B)
  beta1.B <- sample(estimates[,"beta1"],B)
  
  power.bayes <- rep(NA, B)
  
  # Calculate power for each set of sampled values
  for(b in 1:B){
    power.bayes[b] = CalcPower(alpha = alpha, 
                               n = n_village, 
                               m = n_ind, 
                               pii = pii.B[b],
                               tau = tau.B[b],
                               delta = delta.B[b],
                               rho = rho_y.B[b],
                               P0 = exp(beta0.B[b])/(1+exp(beta0.B[b])),
                               P1 = exp(beta0.B[b]+beta1.B[b])/(1+exp(beta0.B[b]+beta1.B[b]))
                               )
  }
  
  # Save power distribution
  filename <- paste("power_", run_name, ".RData", sep ="" )
  save(power.bayes, file = filename)
  
  # Return power distribution
  return(power.bayes)
  
}

