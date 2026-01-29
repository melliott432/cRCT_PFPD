
library(viridis)
library(tidyverse)
library(ggplot2)

colors <- viridis(4)

################################################################################
##### Sample Size Sims #########################################################

setwd("D:/UCSD/Thesis/Bayesian cRCT design/code/sims_sample_size")

n_clusters = c(6, 10, 16, 20)
m <- c(25, 50, 100, 200)

for(j in 1:length(n_clusters)){
  
  power.vec <- c()
  
  power.freq <- c(NA, length = length(m))
  
  for(i in 1:length(m)){
    
    filename = paste("power/power_2026-01-23_", 
                     n_clusters[j], 
                     "_", m[i],
                     ".RData", 
                     sep = "")
    load(filename)
    
    power.vec <- append(power.vec, as.vector(power[[1]]))
    
    power.freq[i] <- power[[2]]
  }
  
  power.df <- data.frame("power" = power.vec,
                         "m" = as.factor(rep(m, each = 1000)))
  
  point.df = data.frame("power" = power.freq, 
                        "m" = m)
  
  m <- factor(m)
  
  ggplot(power.df, aes(x = power, colour = m)) +
    geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75) +
    geom_point(x = point.df$power[1], y = 0, colour = colors[1]) + 
    geom_point(x = point.df$power[2], y = 0, colour = colors[2]) +
    geom_point(x = point.df$power[3], y = 0, colour = colors[3]) +
    geom_point(x = point.df$power[4], y = 0, colour = colors[4]) +
    labs(title = paste("Kernel Density Plot of Power, ", n_clusters[j],  " clusters", sep = ""),
         x = "Power",
         y = "Density")+
    xlim(0,1) + 
    theme_light() +
    scale_color_manual(values = colors)
  ggsave(filename = paste("power_n_", n_clusters[j], ".png", sep =""))
  
}

################################################################################
##### Param Combos Sims ########################################################

setwd("E:/UCSD/Thesis/Bayesian cRCT design/code/sims_param_combos")

n_village = 10
n_ind = 50

pii = seq(0.1, 0.95, by = 0.05)
tau = seq(0.1, 0.95, by = 0.05)
delta = seq(0.1, 0.95, by = 0.05)
sim_params <- expand.grid(pii = pii, tau = tau, delta = delta)
sim_params <- sim_params[1:896,]
power.vec <- c()
power.med <- c()

for(i in 1:nrow(sim_params)){
  
  filename = paste(paste("power/power", 
                         "2026-01-28", 
                         sim_params$pii[i]*100, 
                         sim_params$tau[i]*100, 
                         sim_params$delta[i]*100, 
                         "n_10", 
                         "m_50", 
                         sep = "_"), 
                   "RData", 
                   sep = ".")
  load(filename)
  
  power.vec <- append(power.vec, as.vector(power[[1]]))
  power.med <- append(power.med, median(power[[1]]))
  
}

power.df <- data.frame("power" = power.vec,
                       "pii" = rep(sim_params$pii, each = 1000),
                       "tau" = rep(sim_params$tau, each = 1000),
                       "delta" = rep(sim_params$delta, each = 1000))

power.med <- data.frame("power_median" = power.med,
                        "pii" = sim_params$pii,
                        "tau" = sim_params$tau,
                        "delta" = sim_params$delta)

power.med$delta <- as.character(power.med$delta)
power.med$delta <- as.numeric(power.med$delta)

######### PI vs. TAU #########
######## delta = 0.2 #########

df <- power.med[which(power.med$delta==0.1),]
ggplot(df, aes(pii, tau)) +
  geom_tile(aes(fill = power_median)) +
  scale_fill_viridis(name = "Median Power") +
  ggtitle("Delta = 0.10") +
  xlab("Follow up rate") +
  ylab("Attrition ICC")
ggsave(filename = "plots/power_heatmap_delta10.png")

df <- power.med[which(power.med$delta == 0.15),]
ggplot(df, aes(pii, tau)) +
  geom_tile(aes(fill = power_median)) +
  scale_fill_viridis(name = "Median Power") +
  ggtitle("Delta = 0.15") +
  xlab("Follow up rate") +
  ylab("Attrition ICC")
ggsave(filename = "plots/power_heatmap_delta15.png")

df <- power.med[which(power.med$delta==0.20),]
ggplot(df, aes(pii, tau)) +
  geom_tile(aes(fill = power_median)) +
  scale_fill_viridis(name = "Median Power") +
  ggtitle("Delta = 0.20") +
  xlab("Follow up rate") +
  ylab("Attrition ICC")
ggsave(filename = "plots/power_heatmap_delta20.png")

######### TAU vs. DELTA #########
#########   pi = 0.50   #########

df <- power.med[which(power.med$pii == 0.5),]

ggplot(df, aes(tau, delta)) +
  geom_tile(aes(fill = power_median)) +
  scale_fill_viridis(name = "Median Power") +
  xlab("Attrition ICC") +
  ylab("Treatment effect size")
#ggsave(filename = paste("plots/power_heatmap_pi50_n", n_village, "_m", n_ind, ".png", sep = ""))

######### PI vs. DELTA #########
#########  tau = 0.50  #########

df <- power.med[which(power.med$tau == 0.50 ),]

ggplot(df, aes(pii, delta)) +
  geom_tile(aes(fill = power_median)) +
  scale_fill_viridis(name = "Median Power") +
  xlab("Follow up rate") +
  ylab("Treatment effect size")
#ggsave(filename = paste("plots/power_heatmap_tau50_n", n_village, "_m", n_ind, ".png", sep = ""))



################################################################################
#### Old code from testing
################################################################################

# Power
png(filename = paste(run_name, "_power", ".png", sep = ""))
plot(density(power), main = paste("n = ", n_village, ", m = ", n_ind, sep = ""))
dev.off()

# Power
ggplot(as.data.frame(power), aes(x=power)) + 
  geom_density() + 
  ggtitle(label = paste("n = ", n_village, ", m = ", n_ind, sep = ""))
ggsave(filename = paste(paste(run_name, "power", sep = "_"), "png", sep = "."))

# Power
png(filename = paste("plots/power/", run_name, "_power", ".png", sep = ""))
d = density(power)
d$y = d$y/max(d$y)
plot(d, main = paste("n = ", n_village, ", m = ", n_ind, sep = ""), xlab = "Power")
points(x = freqpwr, y = 0, col = "red", pch = 16)
dev.off()

################################################################################
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