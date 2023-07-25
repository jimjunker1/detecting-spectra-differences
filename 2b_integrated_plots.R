#2b_plots for integrated sims

# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
library(tidybayes) # for nicer looking plots
source("custom_functions.R")
library(cowplot)


# load simulation data ----------------------------------------------------

steep_lambda <- readRDS("data_sim/steep_lambda_dat.rds")
est_lambda <- readRDS("data_sim/est_lambda_dat.rds")
med_lambda <- readRDS("data_sim/med_lambda_dat.rds")
rel_025 <- readRDS("data_sim/relationship_025_dat.rds")
rel_05 <- readRDS("data_sim/relationship_05_dat.rds")
rel_0 <- readRDS("data_sim/relationship_0_dat.rds")
shallow_lambda <- readRDS("data_sim/shallow_lambda_dat.rds")

angle_A <- steep_lambda %>%
  calc_relationship_estimate(.)

angle_B <- med_lambda %>%
  calc_relationship_estimate(.)

angle_C <- shallow_lambda %>%
  calc_relationship_estimate(.)


angles <- bind_rows(steep = angle_A,
                    medium = angle_B,
                    shallow = angle_C,
                    .id = "id")

saveRDS(angles, file = "data_sim/angles.rds")

##### Need to chnage order of legend to be L2n, ELBn, MLE ####


(lambda_window <- angles %>%
  mutate(
    Model = case_when(name == "NAS" ~ "L2n", .default = name),
    Model = factor(Model,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "L2n")),
         sim = factor(id, 
                      levels = c("steep", 
                      "medium", 
                      "shallow"))) %>%
  ggplot(aes(x = estimate, 
             y = Model,
             fill = Model))+
  stat_halfeye(.width = c(0.66, 0.95)) +
  scale_fill_manual(
    values = c("#FF1984",
               "#FF914A",
               "#019AFF"),
    breaks = c("L2n",
               "ELBn", 
               "MLE")) +
  theme_bw() +
  geom_vline(
    aes(xintercept = known_relationship),
    linetype = "dashed") +
  labs(x = "Relationship estimate") +
  facet_wrap(sim~., 
             ncol = 1) +
  NULL)

ggsave(plot = lambda_window,
       filename = "figures/lambda_angle_plot.png",
       units = "in", 
       height = 6,
       width = 6)

rel_data <- bind_rows(beta_0 = rel_0,
                      beta_025 = rel_025,
                      beta_05 = rel_05,
                      .id = "id")

saveRDS(rel_data, file = "data_sim/rel_data.rds")


###### Flip facet-wrap rows, top = 0, bottom = -0.5
rel_data %>%
  mutate(
    Model = case_when(name == "NAS" ~ "L2n", .default = name),
    Model = factor(Model,
                   levels = 
                     c("MLE",
                       "ELBn", 
                       "L2n")),
    # set levels of known relaitonship for facet_wrap() plotting order below
    known_relationship = factor(known_relationship,
                                levels = c("0", "-0.25", "-0.5"))) %>%
  filter(rep <= 100) %>%
  ggplot(aes(x = env_gradient,
             y = estimate, 
             group = rep,
             color = Model)) +
  scale_color_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.05,
              se = FALSE, 
              color = "black") +
  geom_point(alpha = 0.1) +
  facet_wrap(known_relationship~Model)+
  theme_bw() +
  theme(legend.position="none")+
  labs(x = "Hypothetical environemntal gradient")
  NULL

ggsave(filename = "figures/vary_beta_plot.png",
       units = "in", 
       height = 6,
       width = 6)




rel_data_summary = rel_data %>% calc_relationship_estimate(.)

saveRDS(rel_data_summary, file = "data_sim/rel_data_summary.rds")

rel_data_summary %>%
  mutate(
    
    Model = case_when(name == "NAS" ~ "L2n", .default = name),
    Model = factor(Model,
                   levels = 
                     c("MLE",
                       "ELBn", 
                       "L2n"))) %>%
  ggplot(aes(x = estimate, 
             y = Model,
             fill = Model))+
  stat_halfeye(.width = c(0.66, 0.95)) +
  scale_fill_manual(
    values = c("#FF1984",
               "#FF914A",
               "#019AFF"),
    breaks = c("L2n",
               "ELBn", 
               "MLE")) +
  theme_bw() +
  geom_vline(
    aes(xintercept = known_relationship),
    linetype = "dashed") +
  labs(x = "Relationship estimate") +
  facet_wrap(known_relationship~., 
             ncol = 1) +
  NULL

ggsave(filename = "figures/vary_beta_density_plot.png",
       units = "in", 
       height = 6,
       width = 6)

## export one plot as pdf for inkscape practice
rel_data_summary %>%
  filter(known_relationship == -0.5) %>%
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(aes(x = estimate, 
             y = Model,
             fill = Model))+
  stat_halfeye(.width = c(0.66, 0.95)) +
  scale_fill_manual(
    values = c("#FF1984",
               "#FF914A",
               "#019AFF"),
    breaks = c("L2n",
               "ELBn", 
               "MLE")) +
  theme_bw() +
  geom_vline(
    aes(xintercept = known_relationship),
    linetype = "dashed") +
  labs(x = "Relationship estimate") +
  facet_wrap(known_relationship~., 
             ncol = 1) +
  NULL 

# export this figure using the plot window

est_lambda %>%
  mutate(
    Model = case_when(name == "NAS" ~ "L2n", .default = name),
    Model = factor(Model,
                   levels = 
                     c("MLE",
                       "ELBn", 
                       "L2n"))) %>%
  ggplot(
    aes(x = estimate, 
        y = Model,
        fill = Model)) +
  stat_halfeye(.width = c(0.66, 0.95)) +
  scale_fill_manual(
    values = c("#FF1984",
               "#FF914A",
               "#019AFF"),
    breaks = c("L2n",
               "ELBn", 
               "MLE")) +
  theme_bw() +
  geom_vline(aes(xintercept = known_b),
             linetype = "dashed") +
  labs(
    x = "Lambda estimate") +
  facet_wrap(~known_b, scales = "free") +
  scale_y_discrete(breaks=c("MLE",
                            "ELBn", 
                            "NAS"),
                   labels= NULL)

ggsave(filename = "figures/est_lambda_est_b_density.png",
       units = "in", 
       height = 6,
       width = 6)
