############################################################################
#### PLOTS #################################################################

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code/sims_sample_size/plots")

# Power
png(filename = paste(run_name, "_power", ".png", sep = ""))
plot(density(power), main = paste("n = ", n_village, ", m = ", n_ind, sep = ""))
dev.off()

# Power
ggplot(as.data.frame(power), aes(x=power)) + 
  geom_density() + 
  ggtitle(label = paste("n = ", n_village, ", m = ", n_ind, sep = ""))
ggsave(filename = paste(paste(run_name, "power", sep = "_"), "png", sep = "."))

# Follow-up rate
png(filename = paste(paste(run_name, "pii", "posterior", sep = "_"), "png", sep = "."))
coda::densplot(as.mcmc(c(post1[,"pii"], post2[,"pii"], post3[,"pii"], post4[,"pii"])), xlim = c(0,1), show.obs = F, main = paste("n = ", n_village, ", m = ", n_ind, ", Pi Posterior", sep = ""))
points(x = pi.sim, y = 0, col = "red", pch = 16)
dev.off()

# Attrition ICC
png(filename = paste(paste(run_name, "tau", "posterior", sep = "_"), "png", sep = "."))
coda::densplot(as.mcmc(c(post1[,"tau"], post2[,"tau"], post3[,"tau"], post4[,"tau"])), xlim = c((-1/(n_village-1)),1), show.obs = F, main = paste("n = ", n_village, ", m = ", n_ind, ", Tau Posterior", sep = ""))
points(x = tau.sim, y = 0, col = "red", pch = 16)
dev.off()

# Treatment effect
png(filename = paste(paste(run_name, "delta", "posterior", sep = "_"), "png", sep = "."))
coda::densplot(as.mcmc(c(post1[,"delta"], post2[,"delta"], post3[,"delta"], post4[,"delta"])), xlim = c(-1,1), show.obs = F, main = paste("n = ", n_village, ", m = ", n_ind, ", Delta Posterior", sep = ""))
points(x = delta.sim, y = 0, col = "red", pch = 16)
dev.off()

# a 
png(filename = paste(paste(run_name, "a", "posterior", sep = "_"), "png", sep = "."))
coda::densplot(as.mcmc(c(post1[,"a"], post2[,"a"], post3[,"a"], post4[,"a"])), xlim = c(0,3), show.obs = F, main = paste("n = ", n_village, ", m = ", n_ind, ", a Posterior", sep = ""))
points(x = a.sim, y = 0, col = "red", pch = 16)
dev.off()

# b
png(filename = paste(paste(run_name, "b", "posterior", sep = "_"), "png", sep = "."))
coda::densplot(as.mcmc(c(post1[,"b"], post2[,"b"], post3[,"b"], post4[,"b"])), xlim = c(0,3), show.obs = F, main = paste("n = ", n_village, ", m = ", n_ind, ", b Posterior", sep = ""))
points(x = b.sim, y = 0, col = "red", pch = 16)
dev.off()

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code/sims_param_combos")

############################################################################