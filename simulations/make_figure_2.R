
library(grobblR)
library(gridExtra)
library(dplyr)
library(png)
library(viridis)
library(tidyverse)
library(ggplot2)
library(ggpubr)

source("E:/UCSD/Thesis/Bayesian cRCT design/code/function_CalcPower.R", echo = TRUE)

setwd("E:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

################################################################################
# Load in power results and record median and frequentist power

power.df <- data.frame(parms_grid, "med_power" = rep(NA,nrow(parms_grid)), "freq_power" = rep(NA,nrow(parms_grid)))

for(i in 1:nrow(parms_grid)){
    
    run_name = paste(parms_grid[i,"n_village"],parms_grid[i,"n_ind"],100*parms_grid[i,"pii"],100*abs(parms_grid[i,"tau"]),100*abs(parms_grid[i,"delta"]),parms_grid[i,"rep"],sep = "_")
  
    load(file = paste("power_", run_name, ".RData", sep = ""))
    
    power.df$med_power[i] <- median(na.omit(power.dist$power))
    
    power.df$freq_power[i] <- CalcPower(alpha = parms_grid$alpha[i],
                                        n = parms_grid$n_village[i],
                                        m = parms_grid$n_ind[i],
                                        n_arms = parms_grid$n_arm[i],
                                        pii = parms_grid$pii[i],
                                        tau = parms_grid$tau[i],
                                        rho = mean(power.dist$rho),
                                        p1 = parms_grid$p1[i],
                                        p2 = parms_grid$p1[i] + parms_grid$delta[i])
  
}

power.df.1 <- power.df |>
  group_by(pii, tau, delta) |>
  summarise(avg_med_power = mean(med_power),
            avg_freq_power = mean(freq_power))

save(power.df.1, file = "power_results_summary.RData")
################################################################################
## Plot heat maps

load("power_results_summary.RData")

power.df.1 <- power.df.1 |> filter(delta %in% seq(0.10,0.70, by=0.10))

## PI vs. TAU

p.b <- ggplot(power.df.1, aes(pii, tau)) +
  geom_tile(aes(fill = avg_med_power)) +
  scale_fill_viridis(name = "Power",limits=c(0, 1)) +
  labs(title = "") +
  xlab("Follow-up rate") +
  ylab("Attrition ICC") +
  theme(legend.position = "none")

panel.b <- p.b + facet_wrap( ~ delta, nrow = 1, scales = "fixed") + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 90, hjust = 1)) +
  ggtitle("Bayesian Push-forward Power", subtitle = "Treatment Effect")

## PI vs. TAU


p.f<- ggplot(power.df.1, aes(pii, tau)) +
  geom_tile(aes(fill = avg_freq_power)) +
  scale_fill_viridis(name = "Power") +
  labs(title = "") +
  xlab("Follow-up rate") +
  ylab("Attrition ICC") +
  theme(legend.position = "none")

panel.f <- p.f + facet_wrap( ~ delta, nrow = 1,  scales = "fixed") +  
  theme(legend.position = "none", axis.text.x = element_text(angle = 90, hjust=1)) +
  ggtitle("Frequentist Power", subtitle = "Treatment Effect") 


ggarrange(panel.b, panel.f, nrow = 2, ncol = 1,
          common.legend = TRUE, legend = "bottom")
ggsave(filename = "figure2.eps", device=cairo_ps, width = 8, height = 8, units = "in")

setEPS()
postscript(file = "figure2.eps", height = 7, width = 6)
grid.arrange(panel.b, panel.f, nrow = 2, ncol = 1)
dev.off()

?ggarrange()
