# 2 mean and SD estimate plots

library(tidyverse)
source("custom_functions.R")
library(ggpubr)
plb <- read_csv("results/PLB_sim.csv")
plb
ggplot(plb, 
       aes(x = name,
           y = mu_slope, 
           ymin = mu_slope - sd_slope, 
           ymax = mu_slope + sd_slope))+
  geom_pointrange()

ggplot(plb, 
       aes(y = name,
           color = name,
           x = mu_slope, 
           xmin = mu_slope - sd_slope, 
           xmax = mu_slope + sd_slope))+
  geom_pointrange(size = 2) +
  geom_errorbar(inherit.aes = FALSE,
                aes(y = name, 
                    color = name, 
                    xmin = p25, 
                    xmax = p975)) +
  scale_color_viridis_d(option = "plasma") +
  labs(x = "Slope estimate",
       y = "Method") +
  theme_bw() +
  theme(legend.position="none")

plot_mean_slope_estimate <- function(dat){
  ggplot(dat, 
         aes(y = name,
             color = name,
             x = p50, 
             xmin = p25, 
             xmax = p975))+
    geom_pointrange() +
    scale_color_viridis_d(option = "plasma") +
    labs(x = "Slope estimate",
         y = "Method") +
    theme_bw() +
    theme(legend.position="none")
}

(pa <- plot_mean_slope_estimate(plb)+
  geom_vline(xintercept =  -0.5, linetype = "dashed"))
ggsave(file = "figures/plb_slope_estimate.png",
       height = 3,
       width = 5)

pb <- plot_mean_slope_estimate(plb_10)+
  geom_vline(xintercept =  -0.5, linetype = "dashed")
ggsave(file = "figures/plb_10_slope_estimate.png",
       height = 3,
       width = 5)

pc <- plot_mean_slope_estimate(plb_small)+
  geom_vline(xintercept =  -0.1, linetype = "dashed")
ggsave(file = "figures/plb_small_slope_estimate.png",
       height = 3,
       width = 5)

pd <- plot_mean_slope_estimate(plb_static_b) +
  geom_vline(xintercept =  0, linetype = "dashed")
ggsave(file = "figures/plb_static_slope_estimate.png",
       height = 3,
       width = 5)


# making av multi panel plot
ggarrange(pa, pb, pc, pd,
          labels = c("A", "B", "C", "D"))
ggsave("figures/slopes_combined.png",
       width = 5,
       height = 5)

ggplot(plb, 
       aes(y = name,
           x = p50, 
           xmin = p50 - sd_slope, 
           xmax = p50 + sd_slope))+
  geom_pointrange()

ggplot(plb, 
       aes(y = name,
           x = mu_slope, 
           xmin = p25, 
           xmax = p975))+
  geom_pointrange()

plb_static_b <- read_csv("results/plb_static_b.csv")
plb_static_b
ggplot(plb_static_b, 
       aes(x = name,
           y = mu_slope, 
           ymin = mu_slope - sd_slope, 
           ymax = mu_slope + sd_slope))+
  geom_pointrange()

plb_10 <- read_csv("results/PLB_10_sites.csv")
ggplot(plb_10, 
       aes(x = name,
           y = mu_slope, 
           ymin = mu_slope - sd_slope, 
           ymax = mu_slope + sd_slope))+
  geom_pointrange()

plb_small <- read_csv("results/PLB_sim_small.csv")
