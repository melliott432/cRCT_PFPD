
library(R2jags)
library(rjags)
library(gridExtra)
library(truncnorm)
library(tidyverse)
library(grid)
library(ggplot2)
library(beepr)
library(viridis)


setwd("D:/UCSD/Thesis/Bayesian cRCT design")

source("code/function_Estimate.R", echo=TRUE)
source("code/function_GetPowerDistn.R", echo=TRUE)
source("code/function_CalcPower.R", echo=TRUE)

################################################################################
### Format Myanmar data

# Trial data frame columns: 
# "arm_id" ----------------- binary  # 1 = control, 2 = treatment
# "village_id" ------------- integer # 1, 2, 3, 4, 5, 6, ...
# "adherence_final_indic" -- binary  # 0 = non-adherence, 1 = adherence
# "followup_final_indic" --- binary  # 0 = missing, 1 = observed

load("data/data_Bangladesh_pwr.RData")

df <- Bangladesh_pwr.df
df$village_id <- as.numeric(df$village_id)
df <- df %>%
  arrange(village_id)

################################################################################
### Run estimation

run_name = "myanmar"
estimates <- Estimate(run_name = run_name, df = df)
estimates 

b0 <- mean(estimates[,"beta0"])
b1 <- mean(estimates[,"beta1"])

save(estimates, file = "results/application/Bangladesh_estimates.RData")
################################################################################
### Power calculations

# Set operating characteristics and param estimates

alpha = 0.05   # type I error rate
n_ind = 15     # number of individuals per cluster
n_village = 16 # number of clusters total

parms = list("alpha" = alpha,                              
             "m" = n_ind ,                             
             "n" = n_village,
             "n_arms" = 2)         

b0 <- median(estimates[,"beta0"])
b1 <- median(estimates[,"beta1"])

# Calculate Bayesian push-forward power distribution
power <- GetPower(run_name = "Bangladesh", estimates = estimates, parms = parms)

power.bayes <- power$power

power.freq <- CalcPower(alpha = 0.05, n = n_village, m = n_ind, n_arms = 2,
                        pii = median(estimates[,"pii"]), 
                        tau = median(estimates[,"tau"]), 
                        rho = median(estimates[,"rho_y"]), 
                        p1 = exp(b0)/(1+exp(b0)), 
                        p2 = exp(b0+b1)/(1+exp(b0+b1)))

################################################################################

# Save results

save(power.bayes, file = "results/application/Bangladesh_power_bayes.RData")
save(power.freq, file = "results/application/Bangladesh_power_freq.RData")

################################################################################

setwd("D:/UCSD/Thesis/Bayesian cRCT design/results")

estimates <- as.data.frame(estimates)
power.bayes <- as.data.frame(power.bayes)

colors <- viridis(5)

################################################################################
#### Plot predicted power

png("application/Bangladesh_power_figure.png")
ggplot(power.bayes, aes(x = power.bayes)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[4], fill = colors[4]) +
  labs(title = paste("Push-forward Power Distribution for Bangladesh Study", sep = ""),
       subtitle = paste("total number of clusters = ", n_village,", cluster size = ", n_ind, sep = ""),
       x = "Power",
       y = "Density") +
  geom_point(aes(x=power.freq, y=0), color = colors[1], cex = 3) +
  geom_vline(xintercept = median(power.bayes$power.bayes), color = colors[2], lwd = 1, lty = 2) +
  xlim(0,1) +
  theme(legend.position = "right", legend.text = element_text(size = 1, colour ="black"))
dev.off()

################################################################################
#### Plot posterior density and median for parameter estimates

p1 <- ggplot(estimates, aes(x = pii)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[2], fill = colors[2]) +
  labs(x = "Follow-up rate", y = "Density") +
  geom_vline(xintercept = median(estimates$pii),  lty = 2, lwd = 1, color = colors[2]) +
  xlim(0,1) + 
  theme_light()

p2 <- ggplot(estimates, aes(x = tau)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[3], fill = colors[3]) +
  labs(x = "Attrition ICC", y = "Density") +
  geom_vline(xintercept = median(estimates$tau), lty = 2, lwd = 1, color = colors[3]) +
  xlim(0,1) + 
  theme_light()

p3 <- ggplot(estimates, aes(x = delta)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[4], fill = colors[4]) +
  labs(x = "Treatment effect", y = "Density") +
  geom_vline(xintercept = median(estimates$delta), lty = 2, lwd = 1, color = colors[4]) +
  xlim(-1,1) + 
  theme_light()

png("application/Bangladesh_estimates_figure.png", width = 900, height = 300)
grid.arrange(p1, p2, p3,
             nrow = 1, ncol = 3,
             top = "Posterior Distributions")
dev.off()


