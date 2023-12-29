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

# figure 3 ####
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


# figure 4 ####
###### Flip facet-wrap rows, top = 0, bottom = -0.5

ab_data <- data.frame(
  slope = c(0, -0.25, -0.5),
  intercept = -1.5,
  known_relationship = c(0, -0.25, -0.5)
)

ab_data <- ab_data %>%
  expand_grid(method = c("MLE", "ELBn", "L2N"))

rel_data %>%
  filter(rep<=100) %>%
  mutate(
    Model = case_when(name == "NAS" ~ "L2n", .default = name),
    Model = factor(Model,
                   levels = 
                     c("L2n",
                       "ELBn",
                       "MLE"
                       )),
    # set levels of known relaitonship for facet_wrap() plotting order below
    known_relationship = factor(known_relationship,
                                levels = c("0", "-0.25", "-0.5"))) %>%
  filter(rep <= 500) %>%
  ggplot(aes(x = env_gradient,
             y = estimate, 
             group = rep,
             color = Model)) +
  scale_color_manual(
    values = c("#FF1984", "#FF914A", "#019AFF")) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE, 
              color = "black") +
  geom_point(alpha = 0.1) +
  facet_wrap(known_relationship~Model)+
  theme_bw() +
  theme(legend.position="none") +
  labs(x = "Hypothetical environemntal gradient") +
  geom_abline(data = ab_data,
              aes(intercept = intercept, slope = known_relationship),
              color = "red",
              linetype = "dashed",
              linewidth = 1.25,
              alpha = 0.75)+
  NULL

ggsave(filename = "figures/vary_beta_plot.png",
       units = "in", 
       height = 6,
       width = 6)




rel_data_summary = rel_data %>% calc_relationship_estimate(.)

saveRDS(rel_data_summary, file = "data_sim/rel_data_summary.rds")

# Fig 5 ####
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


# Figure 2 #### 
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
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  scale_y_discrete(breaks=c("MLE",
                            "ELBn", 
                            "NAS"),
                   labels= NULL)

ggsave(filename = "figures/est_lambda_est_b_density.png",
       units = "in", 
       height = 6,
       width = 6)

# top panel of figure 2 ####
# add as an SI figure
est_lambda %>%
  mutate(
    Model = case_when(name == "NAS" ~ "L2n", .default = name),
    Model = factor(Model,
                   levels = 
                     c("MLE",
                       "ELBn", 
                       "L2n"))) %>%
  filter(known_b == -2.5 |
           known_b == -2.25 |
           known_b == -2 ) %>%
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
  scale_x_continuous(labels = scales::label_number(accuracy = 0.1)) +
  scale_y_discrete(breaks=c("MLE",
                            "ELBn", 
                            "NAS"),
                   labels= NULL)
ggsave(filename = "figures/SI_fig2_top_panel.png",
       units = "in", 
       height = 6,
       width = 6)


# replicates for SI #### 

angles2 <- angles %>%
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
  filter(known_relationship == -0.5)

# subsample pipeline
angles2 %>%
  group_by(sim) %>%
  filter(rep <= 10) %>%
  mutate(rep_n = 10)

angles_list <- list()

rep_n <- c(10, 50, 100, 200, 250, 500, 750, 1000)
for (i in 1:length(rep_n)){
  j = rep_n[i]
  out <- angles2 %>%
    group_by(sim) %>%
    filter(rep <= j) %>%
    mutate(rep_n = j)
    angles_list[[i]] <- out
}

angles_list %>%
  bind_rows() %>%
  ggplot(aes(x = rep_n,
             y = estimate,
             color = Model))+
  stat_pointinterval(position = "dodge") +
  scale_color_manual(
    values = c("#FF1984",
               "#FF914A",
               "#019AFF"),
    breaks = c("L2n",
               "ELBn", 
               "MLE")) +
  theme_bw() +
  geom_hline(
    aes(yintercept = known_relationship),
    linetype = "dashed") +
  labs(x = "Number of replicates") +
  facet_wrap(sim~., 
             ncol = 1,
             scales = "free_y") +
  NULL

ggsave(filename = "figures/SI_vary_rep_n.png",
       units = "in", 
       height = 6,
       width = 6)
