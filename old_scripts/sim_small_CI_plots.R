# small variation in lambda CI plot

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# simulation data
dat <- readRDS("data_sim/PLB_sim_small_dat.rds")

names(dat)

# estimate the relationship (beta_1_ estimate) across the gradient
relationship_estimate <- dat %>%
  #filter(rep<10) %>%
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

lambda_diff

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


lambda_diff %>%
  pivot_longer(cols = MLE_NAS:NAS_ELB) %>%
  ggplot(aes(y = value,
             x = rep,
             fill = name,
             color = name)) +
  geom_pointrange(aes(ymin = 0, ymax = 0),
                  alpha = 0.3) +
  theme_bw()

dat %>%
  filter(rep %in% seq(0, 1000, by = 10)) %>%
  ggplot(aes(y = estimate - known_b,
             x = rep,
             fill = name,
             color = name)) +
  geom_point() +
  theme_bw() +
  facet_wrap(~known_b)

# proportion of relationships which are significant
relationship_estimate %>%
  filter(p.value < 0.05) %>%
  group_by(name) %>%
  summarize(n() / 1000)



# plot lm_models
# highlight estimates which are significant
dat %>%
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
  select(rep, name, env_gradient,fit, p.value, significant, negative) %>%
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

ggsave("figures/small_lambda_spaghetti.png")


# how far does the estimated relationship deviate from the known relationship?
relationship_estimate %>%
  mutate(diff = known_relationship - estimate) %>%
  group_by(name) %>%
  summarize(mean(diff), 
            sd(diff))

# difference in estimate and known relationship in significant relationships
relationship_estimate %>%
  filter(p.value < 0.05)%>%
  mutate(diff = known_relationship - estimate) %>%
  group_by(name) %>%
  summarize(mean(diff), 
            sd(diff))

# mean estimate of significant relationships
relationship_estimate %>%
  filter(p.value < 0.05)%>%
  group_by(name) %>%
  summarize(mean(estimate), 
            sd(estimate))

# plot
relationship_estimate %>%
  filter(p.value < 0.05) %>%
  ggplot(aes(y = ..scaled..,
             x = estimate, 
             fill = name)) + 
  geom_density(alpha = 0.5,
               adjust = 2) +
  geom_vline(aes(xintercept = known_relationship),
             linetype = "dashed",
             size = 1)+
  theme_bw() +
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Significant Relationships",
       x = "relationship estimate") +
  NULL

# significant positive??
relationship_estimate %>%
  filter(p.value < 0.05, estimate > 0) 

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
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  unnest(c(conf_low, conf_high)) %>%
  select(-data, -lm_mod) %>%
  mutate(ci0 = conf_low <= known_relationship &
           conf_high >= known_relationship) 

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

dat_conf %>%
  filter(p.value < 0.05)%>%
  group_by(name) %>%
  summarize(prop_true = sum(ci0) / n(),
            prop_sig = n()/1000)
