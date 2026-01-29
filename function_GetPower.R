GetPower <- function(df, estimates, parms, run_name){
  
  alpha = parms$alpha
  n_ind = parms$n_ind
  n_village = parms$n_village
  pii = parms$pii
  tau = parms$tau
  delta = parms$delta
  
  ##############################################################################
  
  ## Bayesian Power ##
  
  pr = df |>
    group_by(village_id) |>
    summarize(prop = sum(adherence_base_indic)/length(adherence_base_indic))
  
  var.p = var(pr$prop)
  exp.p = mean(pr$prop)
  rho_y = var.p / (exp.p*(1-exp.p))
  
  p1 = mean(df$adherence_base_indic)
  p2 = delta + p1
  
  B = 1000
  a.B <- sample(estimates[,"a"], B)
  b.B <- sample(estimates[,"b"], B)
  delta.B <- sample(estimates[,"delta"],B)
  tau.B <- 1/(a.B+b.B+1)
  pii.B <- a.B/(a.B+b.B)
  
  G = 1 + (n_ind*pii.B - 1)*rho_y + (1-pii.B)*(1 + (n_ind-1)*tau.B)*rho_y
  
  term = n_village*n_ind*delta.B^2 *(pii.B/G) *(p1*(1-p1) + (delta.B+p1)*(1-(delta.B+p1)))^-1
  
  power.dist <- pnorm(sqrt(term) - qnorm(1-alpha/2))
  
  ##############################################################################
  
  ## Frequentist Power ##
  
  pr = df |>
    group_by(village_id) |>
    summarize(prop = sum(adherence_base_indic)/length(adherence_base_indic))
  
  var.p = var(pr$prop)
  exp.p = mean(pr$prop)
  rho_y = var.p / (exp.p*(1-exp.p))
  
  p1 = mean(df$adherence_base_indic)
  p2 = delta + p1
  
  G = 1 + (n_ind*pii - 1)*rho_y + (1-pii)*(1 + (n_ind-1)*tau)*rho_y
  
  term = n_village*n_ind*delta^2 *(pii/G) *(p1*(1-p1) + (delta+p1)*(1-(delta+p1)))^-1
  
  power.freq <- pnorm(sqrt(term) - qnorm(1-alpha/2))
  
  ##############################################################################
  
  power <- list("power.dist" = power.dist, "power.freq" = power.freq)
  
  filename <- paste("power_", run_name, ".RData", sep ="" )
  save(power, file = filename)
  
  return(power)
  
}