GetPower <- function(run_name, parms, estimates){

  B = 1000
  
  power.dist <- data.frame( "alpha" = rep(parms$alpha,B),
                            "n" = rep(parms$n_village,B),
                            "m" = rep(parms$n_ind,B),
                            "n_arms" = rep(parms$n_arm,B),
                            "pii" = sample(estimates$pii, B, replace = T),
                            "tau" = sample(estimates$tau, B, replace = T),
                            "rho" = sample(estimates$rho_y, B, replace = T),
                            "delta" = sample(estimates$delta, B, replace = T),
                            "beta0" = sample(estimates$beta0, B, replace = T))
  
  # Transform betas to proportions
  power.dist <- power.dist |>
    mutate(p1 = exp(beta0)/(exp(beta0)+1), p2 = p1 + delta) 
  
  # Censor proportions on [0,1]
  power.dist$p1 <- ifelse(power.dist$p1 > 1, 1, power.dist$p1)
  power.dist$p1 <- ifelse(power.dist$p1 < 0, 0, power.dist$p1)
  power.dist$p2 <- ifelse(power.dist$p2 > 1, 1, power.dist$p2)
  power.dist$p2 <- ifelse(power.dist$p2 < 0, 0, power.dist$p2)
  
  # Calculate power for each posterior sample
  power.dist <- power.dist |>
    mutate(power = CalcPower(alpha = alpha, n = n, m = m, n_arms = n_arms, pii = pii, tau = tau, rho = rho, p1 = p1, p2 = p2))

  return(power.dist)
}

