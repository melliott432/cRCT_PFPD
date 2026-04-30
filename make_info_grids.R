
nreps = 5

alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
sd_mult = 0.05 # standard deviation for outcome model

################################################################################
### SAMPLE SIZE COMBINATIONS ###################################################

n_village = c(6, 10, 16, 20)
n_ind = c(25, 50, 100, 200)

pii = 0.42 # follow-up rate
tau = 0.06 # attrition ICC
delta = 0.2 # effect size

ss <- expand.grid("n_village" = n_village, "n_ind" = n_ind, rep = 1:nreps)

siminfo <- data.frame(ss,
                      "n_arm" = n_arm, 
                      "alpha" = alpha,
                      "sd_mult" = sd_mult,
                      "pii" = pii,
                      "tau" = tau,
                      "delta" = delta,
                      "run_name" = paste(ss$n_village, ss$n_ind, ss$rep, sep = "_"),
                      "run_time" = NA)

save(siminfo, file = "D:/UCSD/Thesis/Bayesian cRCT design/results/simulations/sample_size_info.RData")

################################################################################
### PI/TAU/DELTA COMBINATIONS ##################################################

n_village = 10 
n_ind = 20

pii = seq(0.1, 0.95, by = 0.05)
tau = seq(0.1, 0.95, by = 0.05)
delta = seq(0.1, 0.70, by = 0.05)

siminfo <- expand.grid("pii" = pii, "tau" = tau, "delta" = delta, rep = 1:nreps)

siminfo <- data.frame("n_village" = n_village, 
                      "n_ind" = n_ind, 
                      "n_arm" = n_arm, 
                      "alpha" = alpha,
                      "sd_mult" = sd_mult,
                      siminfo,
                      "run_name" = paste(siminfo$pii*100, siminfo$tau*100, siminfo$delta*100, siminfo$rep, sep = "_"),
                      "run_time" = NA)

save(siminfo, file = "D:/UCSD/Thesis/Bayesian cRCT design/results/simulations/param_combos_info.RData")
