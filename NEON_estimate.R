# Empirical example: AMD

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# compare estimates with empirical data
# data is from Pomeranz et al. 2018 Freshwater Biology
neon_dat <- read_csv("neon_data.csv") 

dw_range = range(neon_dat$dw, na.rm = TRUE)

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


neon_result %>%
  group_by(name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ mat.c, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "mat.c") %>%
  select(-data, -lm_mod)%>%
  write_csv("results/NEON_estimates.csv")

neon_result %>%
  group_by(name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ mat.c, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "mat.c") %>%
  select(-data, -lm_mod) %>%
  ggplot(aes(y = estimate, 
             ymin = estimate - std.error, 
             ymax = estimate + std.error,
             x = name,
             color = name)) +
  geom_pointrange()+
  theme_bw() +
  labs(y = "Relationship estimate",
       title = "Change in CSS parameter across Neon")
ggsave("figures/neon_relationship.png")
