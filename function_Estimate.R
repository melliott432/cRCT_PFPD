Estimate <- function(df, run_name) {
  
  # Required packages
  require(R2jags)
  require(rjags)
  
  ##############################################################################
  #### Format data
  
  not_missing = which(is.na(df$adherence_final_indic) == FALSE) # index which observations are not missing
  n_village = length(unique(df$village_id)) # number of villages
  
  # select non-missing outcome observations
  village <- df$village_id[not_missing] # village IDs
  Y <- df$adherence_final_indic[not_missing] # adherence outcome
  X <- df$arm_id[not_missing] - 1 # arm assignment
  
  # calculate cluster follow-up proportions 
  Ri <- df |>
    group_by(village_id) |> # group by village
    summarise(Rprop = sum(followup_final_indic)/length(followup_final_indic)) # calculate proportion
  Ri <- Ri$Rprop
  Ri <- ifelse(Ri>0.999, 0.999, Ri) # set bounds so JAGS doesn't break
  Ri <- ifelse(Ri<0.001, 0.001, Ri)
  R <- df$followup_final_indic # pull individual follow up indicators
  
  # when clusters are unequal sizes
  # create index for individual follow-up variable
  index <- c()
  index[1] <- 1
  for(i in 1:n_village){index[i] <- min(which(df$village_id == i))}
  index[n_village+1] <- (nrow(df)+1)
  
  jagsData <- list("Y" = Y, # outcome
                   "X" = X, # arm
                   "R" = R, # follow-up
                   "index" = index, # cluster sizes
                   "village" = village, # village id
                   "J" = length(Y), # number of outcome observations
                   "N" = length(unique(df$village_id)), # number of villages
                   "pisq" = pi^2) # 3.14^2 (because JAGS doesn't have this)
  
  ############################################################################
  ### Define model for estimation
  
  model_string <- # model as string for JAGS to read
    
    "model {
        
        ### Likelihood ###
        
        for(j in 1:J){ # Individual level
          Y[j] ~ dbern(ilogit(logitp[j]))
          logitp[j] = beta0 + beta1*X[j] + eta[village[j]] + epsilon[j]
          epsilon[j] ~ dnorm(0, prec_eps)
        }
      
        for(n in 1:N){ # Cluster level
            for(m in index[n]:(index[n+1]-1)){
              R[m] ~ dbern(r[n])
            }
          r[n] ~ dbeta(a, b) T(0.001,0.999)
          eta[n] ~  dnorm(0, prec_eta)
        }
      
        ### Priors ###
        
        beta0 ~ dunif(-1.0, 1.0)
        beta1 ~ dunif(-1.0, 1.0)
        a ~ dunif(0, 10.0)
        b ~ dunif(0, 10.0)
        prec_eps ~ dunif(0, 500)
        prec_eta ~ dunif(0, 500)
        
        delta = (exp(beta0+beta1)/(exp(beta0+beta1)+1)) - (exp(beta0)/(exp(beta0)+1))
        tau = 1/(a+b+1)
        pii = a/(a+b)
        rho_y = (1/prec_eta) / ((1/prec_eta) + (pisq/3))
        
      }"
  
  # Compile and initialize JAGS model object
  model <- jags.model(file = textConnection(model_string), data = jagsData, n.chains = 4)
  
  # Run warm up iterations
  update(model, n.iter = 4000)
  
  # Sample from posterior
  posterior_sample <- coda.samples(model, n.iter = 5000,
                                   variable.names = c("pii", "tau", "delta", "a", "b", "rho_y", "beta0", "beta1"))
  
  # Extract posterior samples
  posterior <- rbind(posterior_sample[[1]], posterior_sample[[2]], posterior_sample[[3]], posterior_sample[[4]])
  
  # Save posterior samples
  filename <- paste("posterior_", run_name, ".RData", sep = "")
  save(posterior, file = filename)
  
  return(posterior)
  
}
