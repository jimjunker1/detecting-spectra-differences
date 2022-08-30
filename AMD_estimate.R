# Empirical example: AMD

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# compare estimates with empirical data
# data is from Pomeranz et al. 2018 Freshwater Biology
amd_dat <- read_csv("site_dw_pca.csv") 

dw_range = range(amd_dat$dw, na.rm = TRUE)

# estimate spectra paramaters across AMD gradient using three methods
(amd_result <- amd_dat %>%
  # mutate(date = as.Date(collectDate)) %>%
  group_by(site, pca1) %>%
  #filter(siteID == "CARI") %>%
  # create list-column
  nest() %>%
  mutate(method_compare =
           map(data,
               compare_slopes,
               rsp_var = "dw",
               dw_range = dw_range)) %>%
  ungroup() %>%
  select(-data) %>%
  unnest(cols = method_compare))

# plot parameter estimates (slope or exponent) across AMD gradient
amd_result %>%
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(aes(x = pca1,
           y = estimate,
           color = Model))+
  geom_point() +
  stat_smooth(method = "lm",
              se = FALSE)+
  scale_color_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  theme_bw() +
  labs(title = "Change in exponent across an AMD gradient", 
       x = "AMD gradient (PCA)",
       y = "Slope/exponent estimate")
# save the plot
ggsave("figures/AMD_plot.png")

# estimate the beta_1_ relationship across the gradient
amd_relationship <- amd_result %>%
  group_by(name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ pca1, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "pca1") %>%
  select(-data, -lm_mod) 

# absolute variation across gradient
amd_relationship %>%
  mutate(abs_change = estimate * 
           (max(amd_dat$pca1) - min(amd_dat$pca1))) %>%
  write_csv("results/AMD_estimates.csv")

# plot the mean \pm STD.E of the estimated relationship across the gradient based on method used. 
amd_relationship %>%
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
  scale_color_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  theme_bw() +
  labs(x = "Relationship estimate",
       y = "Method",
       title = "Mean +/- Std. Error Beta across AMD gradient") +
  theme(legend.position="none")

ggsave("figures/amd_relationship.png",
       height = 3,
       width = 5)





# estimate by sample ------------------------------------------------------
# don't think I'm doing this anymore...
# 
# # need to re-export data from AMD project and include "smaple" number
# 
# (amd_result <- amd_dat %>%
#    # mutate(date = as.Date(collectDate)) %>%
#    group_by(site, pca1) %>%
#    #filter(siteID == "CARI") %>%
#    # create list-column
#    nest() %>%
#    mutate(method_compare =
#             map(data,
#                 compare_slopes,
#                 rsp_var = "dw",
#                 dw_range = dw_range)) %>%
#    ungroup() %>%
#    select(-data) %>%
#    unnest(cols = method_compare))