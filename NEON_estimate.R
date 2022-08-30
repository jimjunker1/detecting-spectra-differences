# Empirical example: AMD

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# compare estimates with empirical data
# data is from Pomeranz et al. 2018 Freshwater Biology
neon_dat <- read_csv("neon_data.csv")

dw_range = range(neon_dat$dw, na.rm = TRUE)

# not using the MLEbin method, as in Pomeranz et al. 2022
# just using the "normal" MLE method to compare
# NEON "Bins" individuals into 1 mm widths, so doesn't drastically alter body size estimates

(neon_result <- neon_dat %>%
    # mutate(date = as.Date(collectDate)) %>%
    group_by(mat.c, ID) %>%
    nest() %>%
    mutate(method_compare =
             map(data,
                 compare_slopes,
                 rsp_var = "dw",
                 dw_range = dw_range)) %>%
    ungroup() %>%
    select(-data) %>%
    unnest(cols = method_compare))

ggplot(neon_result,
       aes(x = mat.c,
           y = estimate,
           color = name))+
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE)+
  theme_bw() +
  labs(title = "Change in exponent across NEON", 
       x = "Mean annual air temp",
       y = "Slope/exponent estimate")

ggsave("figures/NEON_plot.png")


neon_relationship <- neon_result %>%
  group_by(name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ mat.c, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "mat.c") %>%
  select(-data, -lm_mod)

neon_relationship %>%
  mutate(abs_change = estimate * 
           (max(neon_dat$mat.c) - min(neon_dat$mat.c))) %>%
  write_csv("results/NEON_estimates.csv")

neon_relationship %>%
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(aes(x = estimate, 
             xmin = estimate - std.error, 
             xmax = estimate + std.error,
             y = Model,
             color = Model)) +
  geom_pointrange()+
  scale_fill_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  theme_bw() +
  labs(y = "Relationship estimate",
       title = "Mean +/- Std. Error Beta across Neon")

ggsave("figures/neon_relationship.png")



# estimate by sample ------------------------------------------------------

# not sure where to go with this want to compare what the estimates are for each sample (i.e. surber sample), and see how that compares to the avergae across samples, and the estimate when combining all of the data for one estimate. 

# need to re-export NEON data from project and include sample number
# estimate CSS parameter for each sample
# compare "total" estimate with individual sample estimates and mean_sample estimates

(neon_sample_result <- neon_dat %>%
   group_by(mat.c, ID, sample) %>%
   mutate(n = n()) %>%
   filter(n >= 1000) %>% # filter out samples that have < 1000 individuals
   nest() %>%
   mutate(method_compare =
            map(data,
                compare_slopes,
                rsp_var = "dw",
                dw_range = dw_range)) %>%
   ungroup() %>%
   select(-data) %>%
   unnest(cols = method_compare))

ggplot(neon_sample_result,
       aes(x = mat.c,
           y = estimate,
           color = name))+
  geom_point(position = position_jitter(widt = 0.2),
             alpha = 0.5) +
  stat_smooth(method = "lm",
              se = FALSE)+
  theme_bw() +
  labs(title = "Change in exponent across NEON", 
       x = "Mean annual air temp",
       y = "Slope/exponent estimate")

neon_sample_result %>%
  group_by(name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ mat.c, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "mat.c") %>%
  select(-data, -lm_mod)

ggplot(neon_sample_result,
       aes(x = mat.c,
           y = estimate))+
  geom_point(position = position_jitter(widt = 0.2),
             alpha = 0.5) +
  stat_summary(aes(group = ID), fun.data = "mean_se", color = "red", size = 0.5) +
  stat_smooth(method = "lm",
              se = FALSE)+
  theme_bw() +
  geom_point(inherit.aes = FALSE,
             data = neon_result, 
             aes(x = mat.c, 
                 y = estimate),
             color = "blue",
             size = 1) +
  facet_wrap(~name, 
             ncol = 1)

neon_sample_result %>%
  filter(ID == 1) %>%
ggplot(
       aes(x = mat.c,
           y = estimate))+
  geom_point(position = position_jitter(),
             alpha = 0.5) +
  stat_summary(aes(group = ID),
               fun.data = "mean_se",
               color = "red",
               size = 0.5) +
  stat_smooth(method = "lm",
              se = FALSE)+
  theme_bw() +
  geom_pointrange(
    inherit.aes = FALSE,
    data = neon_result %>%
      filter(ID == 1),
    aes(x = mat.c,
        y = estimate,
        ymin = minCI, 
        ymax = maxCI),
    color = "dodgerblue") +
  facet_wrap(~name, 
             ncol = 1) +
  labs(title = "NEON Sample Estimates", 
       subtitle = "Black = sample, red = mean sample, blue = combined",
       y = "Slope/exponent estimate") +
  NULL
