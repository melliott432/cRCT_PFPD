
n_reps = 5 # number of replicates to run

alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
sd_mult = 0.05 # standard deviation for outcome model

################################################################################
### SAMPLE SIZE COMBINATIONS ###################################################

n_village = c(12,14,18) # number of clusters total
n_ind = c(25, 50, 100, 200) # cluster size

pii = 0.42 # follow-up rate
tau = 0.06 # attrition ICC
delta = 0.2 # effect size

################################################################################
### PI/TAU/DELTA COMBINATIONS ##################################################

n_village = 10 
n_ind = 20

p1 = 0.3

pii = seq(0.1, 0.95, by = 0.05)
tau = seq(0.1, 0.95, by = 0.05)
delta = seq(0.1, 0.70, by = 0.05)


