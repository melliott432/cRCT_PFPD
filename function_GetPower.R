GetPower <- function(simdata, estimates, run_name){
  
  pr = simdata |>
    group_by(village_id) |>
    summarize(prop = sum(adherence_base_indic)/length(adherence_base_indic))
  
  var.p = var(pr$prop)
  exp.p = mean(pr$prop)
  rho_y = var.p / (exp.p*(1-exp.p))
  
  p1 = mean(simdata$adherence_base_indic)
  p2 = delta + p1
  
  B = 1000
  a.B <- sample(estimates[,"a"], B)
  b.B <- sample(estimates[,"b"], B)
  delta.B <- sample(estimates[,"delta"],B)
  tau.B <- 1/(a.B+b.B+1)
  pii.B <- a.B/(a.B+b.B)
  
  G = 1 + (n_ind*pii.B - 1)*rho_y + (1-pii.B)*(1 + (n_ind-1)*tau.B)*rho_y
  
  term = n_village*n_ind*delta.B^2 *(pii.B/G) *(p1*(1-p1) + (delta.B+p1)*(1-(delta.B+p1)))^-1
  
  power <- pnorm(sqrt(term) - qnorm(1-alpha/2))
  
  filename <- paste("power/power_", run_name, ".RData", sep ="" )
  save(power, file = filename)
  
  return(power)
  
}