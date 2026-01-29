SimEstimate <- function(df, run_name) {
  
  # Format data
  
  n_village = length(unique(df$village_id))
  
  # Adherence outcome
  Y <- na.omit(df$adherence_final_indic)
  
  # Arm assignment
  X <- df$arm_id - 1
  
  # Cluster follow-up proportion 
  Ri <- df |>
    group_by(village_id) |>
    summarise(Rprop = sum(followup_final_indic)/length(followup_final_indic))
  Ri <- Ri$Rprop
  Ri <- ifelse(Ri>0.999, 0.999, Ri)
  Ri <- ifelse(Ri<0.001, 0.001, Ri)
  
  # Individual follow-up 
  R <- matrix(df$followup_final_indic, ncol = n_village)
  
  # Village/cluster ID
  vid <- matrix(df$village_id, ncol = n_village)
  village <- df$village_id
  
  J = length(Y)
  N = ncol(R)
  M = nrow(R)
  
  jagsData <- list("Y" = Y, "X" = X, "R" = R, "village" = village, "J" = J, "N" = N, "M" = M)
  
  ############################################################################
  #### ESTIMATION ############################################################
  
  # Define model for estimation
  model_string <- 
    
    "model {
        
        ### Likelihood ###
        
        for(j in 1:J){ # Individual level
          Y[j] ~ dbern(ilogit(logitp[j]))
          logitp[j] = beta*X[j] + eta[village[j]] + epsilon[j]
          epsilon[j] ~ dnorm(0, prec_eps)
        }
      
        for(n in 1:N){ # Cluster level
          for(m in 1:M){
            R[m,n] ~ dbern(r[n])
          }
          r[n] ~ dbeta(a, b) T(0.001,0.999)
          eta[n] ~  dnorm(0, prec_eta)
        }
      
        ### Priors ###
        
        beta ~ dunif(-1.0, 1.0)
        a ~ dunif(0, 10.0)
        b ~ dunif(0, 10.0)
        prec_eps ~ dunif(0, 500)
        prec_eta ~ dunif(0, 500)
        
        delta <- (exp(beta)/(exp(beta)+1)) - 0.50
        tau = 1/(a+b+1)
        pii = a/(a+b)
        
      }"
  
  # Compile and initialize JAGS model object
  model <- jags.model(file = textConnection(model_string), 
                      data = jagsData,
                      n.chains = 4)
  
  # Warm up iterations
  update(model, n.iter = 4000)
  
  # Sampling from posterior
  posterior_sample <- coda.samples(model,
                                   variable.names = c("pii", "tau", "delta", "a", "b"), 
                                   n.iter = 5000)
  
  # Extract posterior samples
  post1 <- posterior_sample[[1]]
  post2 <- posterior_sample[[2]]
  post3 <- posterior_sample[[3]]
  post4 <- posterior_sample[[4]]
  posterior <- rbind(post1, post2, post3, post4)
  
  filename <- paste("posterior_", run_name, ".RData", sep = "")
  save(posterior, file = filename)
  
  return(posterior)
  
}