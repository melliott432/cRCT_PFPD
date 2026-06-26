
library(viridis)
library(ggplot2)
library(tidyverse)
library(latex2exp)

n_reps = 5 # number of replicates to run
alpha = 0.05 # significance level
n_arm = 2 # number of experimental arms
sd_mult = 0.05 # standard deviation for outcome model
p1 = 0.3
n_village = c(10,12,14,18) # number of clusters total
n_ind = c(25, 50, 100, 200) # cluster size
pii = 0.42 # follow-up rate
tau = c(0.01, 0.06, 0.1, 0.3) # attrition ICC
delta = 0.2 # effect size

parms_grid <- expand.grid("rep" = 1:n_reps, "alpha" = alpha, "n_arm" = n_arm, 
                          "sd_mult" = sd_mult, "p1" = p1, 
                          "n_village" = n_village, "n_ind" = n_ind, 
                          "pii" = pii, "tau" =  tau, "delta" = delta)

power.df <- data.frame(NULL)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

for(i in 1:nrow(parms_grid)){
  
  run_name = paste(parms_grid[i,"n_village"],parms_grid[i,"n_ind"],100*parms_grid[i,"pii"],100*abs(parms_grid[i,"tau"]),100*abs(parms_grid[i,"delta"]),parms_grid[i,"rep"],sep = "_")
  
  load(file=paste("power_",run_name,".RData", sep = ""))
  
  power.df1 <- data.frame("n" = rep(parms_grid[i,"n_village"]), 
                          "m" = parms_grid[i,"n_ind"],
                          "tau" = parms_grid[i,"tau"],
                          "power" = power.dist$power)
  
  power.df <- rbind(power.df, power.df1)
}


custom_titles <- c(
  "10" = "n = 10",
  "12" = "n = 12",
  "14" = "n = 14",
  "16" = "n = 18",
  "0.01" = TeX(r'($\tau = 0.01$)'),
  "0.06" = TeX(r'($\tau = 0.06$)'),
  "0.1" = TeX(r'($\tau = 0.1$)'),
  "0.3" = TeX(r'($\tau = 0.3$)')
)

ggplot(power.df, aes(x = power, colour = as.factor(m)
                     ,linetype = as.factor(m), fill = as.factor(m)
                     )) +
  facet_wrap(~ n + tau, labeller = labeller(
    n = ~ paste("n =", .),
    tau = ~ paste("ICC =", .),
    .multi_line = FALSE
  )) +
  geom_density(mapping = aes(y = after_stat(scaled)), lwd = 0.5, alpha = 0.2,
               key_glyph = "path") +
  labs(title = "Push-forward Power Distributions",
       x = "Power",
       y = "Density",
       color = "Cluster Size",
       fill = "Cluster Size",
       linetype = "Cluster Size") +
  xlim(0, 1) +
  scale_color_manual(values = viridis(7)) +
  scale_linetype_manual(values = c("twodash",  "dashed", "longdash", "solid")) +
  scale_fill_manual(values = viridis(6)) +
  guides(color = guide_legend(override.aes = list(fill = NA, alpha = NA))) +
  theme(legend.position = "right", axis.text.x = element_text(angle = 90, hjust=1))

ggsave(filename = "figure1.eps", device=cairo_ps, width = 8, height = 8, units = "in")

