# Eqn 13 from Taljaard, Donner, and Klar 2007

CalcPower <- function(alpha, n, m, pii, tau, delta, rho, P0, P1){
  
  G <- 1 + (((m*pii) - 1)*rho) + (1-pii)*(1 + (m-1)*tau)*rho # correction factor
  
  A <- pii*(delta^2)*n*m # numerator
  
  B <- (P1*(1-P1)) + (P0*(1-P0)) # denominator
  
  za2 <- qnorm(alpha/2) # alpha/2 z value
  
  zb <- sqrt(A/(B*2*G)) - za2 # calculate beta z value
    
  power <- 1 - pnorm(zb) # power = 1 - beta

  return(power)
  
}
