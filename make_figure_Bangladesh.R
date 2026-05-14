################################################################################
### Make plots

library(viridis)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/results")

load("D:/UCSD/Thesis/Bayesian cRCT design/results/application/Bangladesh_power_freq.RData")
load("D:/UCSD/Thesis/Bayesian cRCT design/results/application/Bangladesh_power_bayes.RData")
load("D:/UCSD/Thesis/Bayesian cRCT design/results/application/Bangladesh_estimates.RData")

estimates <- as.data.frame(estimates)
power.bayes <- as.data.frame(power.bayes)

colors <- viridis(5)

################################################################################
#### Plot predicted power

png("application/Bangladesh_power_figure.png")
ggplot(power.bayes, aes(x = power.bayes)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[4], fill = colors[4]) +
  labs(title = paste("Density Plot of Power for Bangladesh Study", sep = ""),
       subtitle = "16 clusters of size 15",
       x = "Power",
       y = "Density") +
  geom_point(aes(x=power.freq, y=0), color = colors[1], cex = 3) +
  geom_vline(xintercept = median(power.bayes$power.bayes), color = colors[2], lwd = 1, lty = 2) +
  xlim(0,1) + 
  theme_light()
dev.off()

################################################################################
#### Plot posterior density and median for parameter estimates

p1 <- ggplot(estimates, aes(x = pii)) +
        geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[2], fill = colors[2]) +
        labs(x = "Follow-up rate", y = "Density") +
        geom_point(data = estimates, mapping = aes(x = median(pii), y = 0), pch = 4, cex = 2, stroke = 2, color = colors[2]) +
        xlim(0,1) + 
        theme_light()

p2 <- ggplot(estimates, aes(x = tau)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[3], fill = colors[3]) +
  labs(x = "Attrition ICC", y = "Density") +
  geom_point(data = estimates, mapping = aes(x = median(tau), y = 0), pch = 4, cex = 2, stroke = 2, color = colors[3]) +
  xlim(0,1) + 
  theme_light()

p3 <- ggplot(estimates, aes(x = delta)) +
  geom_density(mapping = aes(y = after_stat(scaled)), alpha = 0.75, linewidth = 0.75, color = colors[4], fill = colors[4]) +
  labs(x = "Treatment effect", y = "Density") +
  geom_point(data = estimates, mapping = aes(x = median(delta), y = 0), pch = 4, cex = 2, stroke = 2, color = colors[4]) +
  xlim(-1,1) + 
  theme_light()

png("application/Bangladesh_estimates_figure.png", width = 400, height = 400)
grid.arrange(p1, p2, p3,
             nrow = 1, ncol = 3,
             top = "Posterior Distributions")
dev.off()

 