# pairwise estimates of lambda

dat <- readRDS("data_sim/PLB_shallow_lambda_dat.rds")

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
  summarize(n() / 1000)

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
  summarize(prop_true = sum(ci0))
# NAS 
# 748 do not have -0.2 in CI
dat_conf %>%
  ungroup() %>%
  filter(ci0 == FALSE, name == "NAS", conf_low < -0.199) %>%
  count()
# 25 have bound at -0.199
# (748 + 25) / 100 = 77.3%


dat_conf %>%
  ungroup() %>%
  filter(ci0 == FALSE, name == "NAS", conf_low < -0.198) %>%
  count()

dat_conf %>%
  group_by(name) %>%
  mutate(conf_width = conf_high - conf_low) %>%
  summarize(mean_width = mean(conf_width),
            sd_width = sd(conf_width))
