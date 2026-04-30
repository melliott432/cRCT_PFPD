
library(viridis)
library(tidyverse)
library(ggplot2)
library(gridExtra)
library(ggpubr)

setwd("D:/UCSD/Thesis/Bayesian cRCT design/results/simulations")

load("param_combos_power_summary.RData")

################################################################################
#### Figure X: Bayes vs. Freq Power Heatmaps ###################################
################################################################################

#### Bayes #####################################################################

df <- power.summary |>
  filter(delta < 0.70)

pii.df <- df |>
  filter(pii == c(0.10, 0.25, 0.50, 0.75, 0.90))

tau.df <- df |>
  filter(tau == c(0.10, 0.25, 0.50, 0.75, 0.90))

delta.df <- df |>
  filter(delta == c(0.10, 0.25, 0.40, 0.55, 0.70))


## PI vs. TAU

p11 <- ggplot(delta.df, aes(pii, tau)) +
        geom_tile(aes(fill = avg_median_power)) +
        scale_fill_viridis(name = "Median Power") +
        labs(title = "") +
        xlab("Follow-up rate") +
        ylab("Attrition ICC") +
        theme(legend.position = "none")

panel1 <- p11 + facet_wrap( ~ delta, nrow = 1, scales = "fixed") + 
                theme(legend.position = "right") +
                ggtitle("Pi vs. Tau (Bayes)")



## TAU vs. DELTA

p21 <- ggplot(pii.df, aes(tau, delta)) +
  geom_tile(aes(fill = avg_median_power)) +
  scale_fill_viridis(name = "Median Power") +
  labs(title = "") +
  xlab("Attrition ICC") +
  ylab("Treatment effect") +
  theme(legend.position = "none")

panel2 <- p21 + facet_wrap( ~ pii, nrow = 1, scales = "fixed") + 
  theme(legend.position = "right") +
  ggtitle("Tau vs. Delta (Bayes)")



## PI vs. DELTA

p31 <- ggplot(tau.df, aes(pii, delta)) +
  geom_tile(aes(fill = avg_median_power)) +
  scale_fill_viridis(name = "Median Power") +
  labs(title = "") +
  xlab("Follow-up rate") +
  ylab("Treatment effect") +
  theme(legend.position = "none")

panel3 <- p31 + facet_wrap( ~ tau, nrow = 1,  scales = "fixed") + 
  theme(legend.position = "right") +
  ggtitle("Pi vs. Delta (Bayes)")


#### Freq ######################################################################

## PI vs. TAU


p41 <- ggplot(delta.df, aes(pii, tau)) +
  geom_tile(aes(fill = freq_power)) +
  scale_fill_viridis(name = "Frequentist Power") +
  labs(title = "") +
  xlab("Follow-up rate") +
  ylab("Attrition ICC") +
  theme(legend.position = "none")

panel4 <- p41 + facet_wrap( ~ delta, nrow = 1,  scales = "fixed") +  
  theme(legend.position = "right") +
  ggtitle("Pi vs. Tau (Freq)")


## TAU vs. DELTA


p51 <- ggplot(pii.df, aes(tau, delta)) +
  geom_tile(aes(fill = freq_power)) +
  scale_fill_viridis(name = "Frequentist Power") +
  labs(title = "") +
  xlab("Attrition ICC") +
  ylab("Treatment effect") +
  theme(legend.position = "none")

panel5 <- p51 + facet_wrap( ~ pii, nrow = 1,  scales = "fixed") + 
  theme(legend.position = "right") +
  ggtitle("Tau vs. Delta (Freq)")

## PI vs. DELTA

p61 <- ggplot(tau.df, aes(pii, delta)) +
  geom_tile(aes(fill = freq_power)) +
  scale_fill_viridis(name = "Frequentist Power") +
  labs(title = "") +
  xlab("Follow-up rate") +
  ylab("Treatment effect") +
  theme(legend.position = "none")

panel6 <- p61 + facet_wrap( ~ tau, nrow = 1,  scales = "fixed") +  
  theme(legend.position = "right") +
  ggtitle("Pi vs. Delta (Freq)")


################################################################################
## Arrange panels and save image

png("param_combos_figure.png", width = 2000, height = 1500)
grid.arrange(panel1, panel4, panel2, panel5, panel3, panel6,
             nrow = 3, ncol = 2)
dev.off()


