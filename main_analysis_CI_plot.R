# main analysis CI plot

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# main simulation data
dat <- readRDS("data_sim/PLB_sim_dat.rds")

names(dat)

# estimate the relationship (beta_1_ estimate) across the gradient
relationship_estimate <- dat %>%
  group_by(rep, name, known_relationship) %>%
  nest() %>%
  mutate(lm_mod =
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  select(-data, -lm_mod, -statistic)

# what is the pairwise difference in spectra parameter estimates?
lambda_diff <- dat %>%
  select(env_gradient, rep, name, estimate) %>%
  pivot_wider(names_from = name, values_from = estimate) %>%
  mutate(MLE_NAS = MLE - NAS,
         MLE_ELB = MLE - ELBn, 
         NAS_ELB = NAS - ELBn)

lambda_diff %>%
  pivot_longer(cols = MLE_NAS:NAS_ELB) %>%
  ggplot(aes(y = ..scaled..,
             x = value,
             fill = name)) +
  geom_density(alpha = 0.55,
               adjust = 2) +
  geom_vline(aes(xintercept = 0),
             linetype = "dashed",
             color = "black",
             size = 1) +
  theme_bw() +
  scale_fill_viridis_d(option = "A", direction = -1)


# proportion of relationships which are significant
relationship_estimate %>%
  filter(p.value < 0.05) %>%
  group_by(name) %>%
  summarize(n())

# how far does the estimated relationship deviate from the known relationship?
relationship_estimate %>%
  mutate(diff = known_relationship - estimate) %>%
  group_by(name) %>%
  summarize(mean(diff), 
            sd(diff))


dat_conf <- dat %>%
  #filter(rep <= 10) %>%
  group_by(rep, name, known_relationship) %>%
  nest() %>%
  mutate(lm_mod =
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(conf_low = map(lm_mod, 
                        ~confint(., parm = 2)[1]),
         conf_high = map(lm_mod, 
                         ~confint(., parm = 2)[2])) %>%
  unnest(c(conf_low, conf_high)) %>%
  select(-data, -lm_mod) %>%
  mutate(ci0 = conf_low <=-0.5 & conf_high >=-0.5) 

dat_conf %>%
  filter(ci0 == FALSE) %>%
  ggplot(aes(x = rep, ymin = conf_low, ymax = conf_high, color = ci0)) +
  geom_linerange() +
  facet_wrap(~name) +
  theme_bw()

dat_conf %>%
  group_by(name) %>%
  summarize(prop_true = sum(ci0) / n())

dat_conf %>%
  group_by(name) %>%
  mutate(conf_width = conf_high - conf_low) %>%
  summarize(mean_width = mean(conf_width),
            sd_width = sd(conf_width))
