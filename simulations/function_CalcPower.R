CalcPower <- function(alpha, n, m, n_arms, pii, tau, rho, p1, p2){
  
  z.a2 = qnorm(alpha/2)
  
  k = n/n_arms

  G <- 1 + (((m*pii)-1)*rho)+(1-pii)*(1+(m-1)*tau)*rho
  
  z.b <- sqrt(k*m*pii*(1/G)*((p1-p2)^2)*(((p2*(1-p2))+(p1*(1-p1)))^(-1))) + z.a2 
  
  power <- pnorm(z.b)
  
  return(power)
}
