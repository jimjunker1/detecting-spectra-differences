# no relationship CI plot

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# no relationship data
dat <- readRDS("data_sim/PLB_static_b_dat.rds")

names(dat)

# calculate beta coefficients for all reps
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

# how many show a "significant" result?
relationship_estimate %>%
  filter(p.value < 0.05) %>%
  group_by(name) %>%
  summarize(n()/1000)

spaghetti_dat <- dat %>%
  #filter(rep < 600) %>%
  group_by(rep, name, known_relationship) %>%
  nest() %>%
  mutate(lm_mod =
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(fit = map(lm_mod,
                   ~fitted(object = .))) %>%
  unnest(c(fit, data)) %>%
  select(rep, name, env_gradient, fit, lm_mod) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  mutate(significant = p.value < 0.05,
         negative = estimate < 0) %>%
  select(-lm_mod, -term, -std.error, -statistic) %>%
  ungroup() %>%
  select(rep, name, env_gradient,fit, p.value, significant, negative) 

spaghetti_dat %>%
  ggplot(aes(y = fit,
             x = env_gradient,
             group = rep, 
             color = negative)) +
  geom_point() +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.25,
              se = FALSE,
              size = 1) +
  facet_wrap(name~significant,
             labeller = label_both,
             ncol = 2) +
  scale_color_viridis_d(option = "plasma",
                        begin = 0,
                        end = 0.6) +
  theme_bw() +
  NULL

ggsave("figures/static_b_lambda_spaghetti.png")

spaghetti_dat %>%
  filter(significant == TRUE) %>%
  group_by(name) %>%
  summarize(sum(negative) / n())

# what is the difference in estimates and known relationship?
relationship_estimate %>%
  mutate(diff = known_relationship - estimate) %>%
  group_by(name) %>%
  summarize(mean(diff), 
            sd(diff))

# confidence intervals of beta estimates
dat_conf <- dat %>%
  group_by(rep, name, known_relationship) %>%
  nest() %>%
  mutate(lm_mod =
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(conf_low = map(lm_mod, 
                    ~confint(., parm = 2)[1]),
         conf_high = map(lm_mod, 
                        ~confint(., parm = 2)[2])) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  unnest(c(conf_low, conf_high)) %>%
  select(-data, -lm_mod) %>%
  mutate(ci0 = conf_low <=0 & conf_high >=0) 

# plot of CI's color-coded if contain tru value or not
# This could probably be improved...
dat_conf %>%
  #filter(ci0 == FALSE) %>%
  ggplot(aes(x = estimate, ymin = conf_low, ymax = conf_high, color = ci0)) +
  geom_linerange() +
  facet_wrap(~name) +
  theme_bw()

# alternate format
dat_conf %>%
  filter(ci0 == FALSE) %>%
  ggplot(aes(x = rep, ymin = conf_low, ymax = conf_high, color = ci0)) +
  geom_linerange() +
  facet_wrap(~name) +
  theme_bw()


# alternate format
dat_conf %>%
  filter(ci0 == TRUE) %>%
  ggplot(aes(x = rep, ymin = conf_low, ymax = conf_high, color = ci0)) +
  geom_linerange() +
  facet_wrap(~name) +
  theme_bw()

# proportion of reps with known relationship in CI
dat_conf %>%
  group_by(name) %>%
  summarize(prop_true = sum(ci0) / n())

# width of CI's
dat_conf %>%
  group_by(name) %>%
  mutate(conf_width = conf_high - conf_low) %>%
  summarize(mean_width = mean(conf_width),
            sd_width = sd(conf_width))

# ~ 95% of estimates have known value in CI regardless of method
# However, mean and SD CIs are ~ 3x wider in binned methods
