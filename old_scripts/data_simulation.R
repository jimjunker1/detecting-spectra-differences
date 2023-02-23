# script 1
# simulate data sets
# analyze using different methods
# Compare results


# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
library(VGAM) # truncated pareto distribution
source("custom_functions.R")


# Figuring out simulation ----------------------------------------------------------------
# 
# 
# n <-  1000
# b <- c(-1.5, -2, -2.5)
# x <- c(-1, 0, 1)
# m_15 <- rPLB(n = n, b = b[1])
# m_20 <- rPLB(n = n, b = b[2])
# m_25 <- rPLB(n = n, b = b[3])
# m_sample <- c(m_15, m_20, m_25)
# 
# df <- tibble(known_b = 
#                rep(b,
#                    each = n),
#              x = rep(x, each = n),
#              n = n,
#              name = 
#                rep(c("m_invert_1.5",
#                      "m_invert_2.0",
#                      "m_invert_2.5"),
#                    each = n))
# df$m <- m_sample
# 
# # just MLE estimate
# # mle_lambda <- df %>%
# #   group_by(known_b) %>%
# #   nest() %>%
# #   mutate(lambda = map(data,
# #                       MLE_tidy,
# #                       "m")) %>%
# #   unnest(cols = lambda) %>%
# #   select(-data) %>%
# #   ungroup()
# 
# # compare slopes method from Pomeranz et al. 2022
dw_range = range(df$m, na.rm = TRUE)
# 
# (test_method <- df %>%
#   # mutate(date = as.Date(collectDate)) %>%
#   group_by(known_b, x, rep) %>%
#   #filter(siteID == "CARI") %>%
#   # create list-column
#   nest() %>% 
#   mutate(method_compare = 
#            map(data,
#                compare_slopes,
#                rsp_var = "m",
#                dw_range = dw_range)) %>%
#   ungroup() %>%
#   select(-data) %>%
#   unnest(cols = method_compare))
# 


# scaled up simulation ----------------------------------------------------

# Scale up simulation ####

n <-  1000
b <- c(-1.5, -2, -2.5)
x <- c(-1, 0, 1)
rep <- 1000
sim_out <- list()
#df <- tibble()

for(i in 1:rep){
  df <- tibble(
    known_b = rep(b, each = n),
    x = rep(x, each = n),
    n = n)
  x_15 <- rPLB(n = n, b = b[1])
  x_20 <- rPLB(n = n, b = b[2])
  x_25 <- rPLB(n = n, b = b[3])
  x_sample <- c(x_15, x_20, x_25)
  df$m <- x_sample
  df$rep <- i
  sim_out[[i]] <- df
}

sim_df <- bind_rows(sim_out)
# this takes a while... surely could optimize code to run faster
(test_method <- sim_df %>%
    # mutate(date = as.Date(collectDate)) %>%
    group_by(rep, known_b, x) %>%
    #filter(siteID == "CARI") %>%
    # create list-column
    nest() %>% 
    mutate(method_compare = 
             map(data,
                 compare_slopes,
                 rsp_var = "m",
                 dw_range = dw_range)) %>%
    ungroup() %>%
    select(-data) %>%
    unnest(cols = method_compare))

# plot each regression "rep"
ggplot(test_method,
       aes(x = x,
           y = estimate, 
           group = rep,
           color = rep)) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE)+
  geom_point(#color = "red"
    ) +
  facet_wrap(~name)+
  theme_bw()

# distribution of slope estimates?
test_method %>%
  group_by(rep, name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ x, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "x") %>%
  group_by(name) %>%
  summarize(mu_slope = mean(estimate, na.rm = TRUE),
            sd_slope = sd(estimate),
            p25 = quantile(estimate, probs = 0.25),
            p975 = quantile(estimate, probs = 0.975))




# simulation function -----------------------------------------------------

# PLB ####
rep = 100
PLB_sim <- sim_result(b = c(-1.5, -2, -2.5),
           env_gradient = c(-1, 0, 1),
           rep = rep, 
           m_lower = 0.0026,
           m_upper = 1.2e4,
           distribution = "PLB")


# plot each regression "rep"
ggplot(PLB_sim,
       aes(x = env_gradient,
           y = estimate, 
           group = rep,
           color = rep)) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE)+
  geom_point() +
  facet_wrap(~name)+
  labs(title = "Bouned Power Law",
       subtitle = "m_range = (0.0026, 12,000) \n x = (-1, 1); b = (-2.5, -1.5)") +
  theme_bw()

ggsave("figures/plb_main.png")

# distribution of slope estimates?
PLB_sim %>%
  group_by(rep, name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  group_by(name) %>%
  summarize(mu_slope = mean(estimate, na.rm = TRUE),
            sd_slope = sd(estimate),
            p25 = quantile(estimate, probs = 0.25),
            p975 = quantile(estimate, probs = 0.975))

PLB_sim %>%
  ggplot(aes(x = estimate, 
             fill = name)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  facet_wrap(~known_b,
             scales = "free") +
  geom_vline(aes(xintercept = known_b),
             size = 1,
             alpha = 0.75,
             color = "black")+
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Bounded Power Law",
       x = "slope estimate") +
  NULL
ggsave("figures/plb_b_est_dist.png")

# truncated pareto ####
tp_sim <- sim_result(b = c(0.1, 0.6, 1.1),
                      env_gradient = c(-1, 0, 1),
                      rep = rep, 
                      m_lower = 0.0026,
                      m_upper = 1.2e4,
                      distribution = "tpareto")


# plot each regression "rep"
ggplot(tp_sim,
       aes(x = env_gradient,
           y = estimate, 
           group = rep)) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE)+
  geom_point(#color = "red"
  ) +
  facet_wrap(~name)+
  labs(title = "Truncated Pareto",
       subtitle = "m_min = 0.0026 \nm_max = 12,000") +
  theme_bw()

# distribution of slope estimates?
tp_sim %>%
  group_by(rep, name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  group_by(name) %>%
  summarize(mu_slope = mean(estimate, na.rm = TRUE),
            sd_slope = sd(estimate),
            p25 = quantile(estimate, probs = 0.25),
            p975 = quantile(estimate, probs = 0.975))
# smaller slope variation -------------------------------------------------
PLB_sim <- sim_result(b = c(-1.9, -2, -2.1),
                      env_gradient = c(-1, 0, 1),
                      rep = rep, 
                      m_lower = 0.0026,
                      m_upper = 1.2e4,
                      distribution = "PLB")


# plot each regression "rep"
ggplot(PLB_sim,
       aes(x = env_gradient,
           y = estimate, 
           group = rep)) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE)+
  geom_point(#color = "red"
  ) +
  facet_wrap(~name)+
  labs(title = "Bouned Power Law: b: -1.9 to -2.1",
       subtitle = "m_min = 0.0026 \nm_max = 12,000") +
  theme_bw()

# distribution of slope estimates?
PLB_sim %>%
  group_by(rep, name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  group_by(name) %>%
  summarize(mu_slope = mean(estimate, na.rm = TRUE),
            sd_slope = sd(estimate),
            p25 = quantile(estimate, probs = 0.25),
            p975 = quantile(estimate, probs = 0.975))

PLB_sim %>%
  ggplot(aes(x = estimate, 
             fill = name)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  facet_wrap(~known_b, 
             scales = "free_x") +
  geom_vline(aes(xintercept = known_b),
             size = 1,
             alpha = 0.75,
             color = "black")+
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Bounded Power Law") +
  NULL

# truncated pareto ####
tp_sim <- sim_result(b = c(0.4, 0.5, .6),
                     env_gradient = c(-1, 0, 1),
                     rep = rep, 
                     m_lower = 0.0026,
                     m_upper = 1.2e4,
                     distribution = "tpareto")


# plot each regression "rep"
ggplot(tp_sim,
       aes(x = env_gradient,
           y = estimate, 
           group = rep)) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE)+
  geom_point(#color = "red"
  ) +
  facet_wrap(~name)+
  labs(title = "Truncated Pareto, shape: 0.4 to 0.6",
       subtitle = "m_min = 0.0026 \nm_max = 12,000") +
  theme_bw()

# distribution of slope estimates?
tp_sim %>%
  group_by(rep, name) %>%
  nest() %>%
  mutate(lm_mod = 
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  group_by(name) %>%
  summarize(mu_slope = mean(estimate, na.rm = TRUE),
            sd_slope = sd(estimate),
            p25 = quantile(estimate, probs = 0.25),
            p975 = quantile(estimate, probs = 0.975))

tp_sim %>%
  ggplot(aes(x = estimate, 
             fill = name)) +
  geom_density(alpha = 0.5) +
  theme_bw() +
  facet_wrap(~known_b, 
             scales = "free_x") +
  geom_vline(
    aes(
      xintercept = (known_b*-1) -1),
             size = 1,
             alpha = 0.75,
             color = "black")+
  scale_fill_viridis_d(option = "plasma") +
  labs(title = "Truncated Pareto") +
  NULL

# other analyses - beta ---------------------------------------------------


test_method %>%
  mutate(diff_est = known_b - estimate) %>%
  group_by(name, known_b) %>%
  summarize(mean_diff = mean(diff_est))

test_method <- test_method %>%
  #filter(rep < 2) %>%
  mutate(in_CI = minCI < known_b & maxCI > known_b)

ggplot(test_method,
       aes(y = estimate, 
           ymin = minCI, 
           ymax = maxCI,
           x = rep)) +
  geom_pointrange() +
  geom_point(inherit.aes = FALSE,
             aes(y = known_b,
                 x = rep), 
             color = "red",
             shape = 17)+
  facet_wrap(~name)

ggplot(test_method,
       aes(y = estimate, 
           ymin = minCI, 
           ymax = maxCI,
           x = rep,
           color = name)) +
  geom_pointrange() +
  geom_point(inherit.aes = FALSE,
             aes(y = known_b,
                 x = rep), 
             color = "red",
             shape = 17)+
  facet_wrap(name~known_b,
             scales = "free")

test_method %>%
  ggplot(aes(y = estimate, 
             ymin = minCI, 
             ymax = maxCI,
             x = rep,
             color = in_CI)) +
  geom_pointrange() +
  facet_wrap(~name)

test_method %>% group_by(name) %>% summarize(mean(in_CI))

test_method %>%
  ggplot(aes(x = known_b,
             y = estimate, 
             color = name)) +
  geom_point() +
  stat_smooth(method = "lm")

coef(lm(estimate ~ known_b, data = test_method %>%
          filter(name == "ELBn")))
coef(lm(estimate ~ known_b, data = test_method %>%
          filter(name == "MLE")))
coef(lm(estimate ~ known_b, data = test_method %>%
          filter(name == "NAS")))


# smaller difference between slopes
n <-  1000
b <- c(-1.9, -2, -2.1)
rep <- 100
sim_out <- list()
#df <- tibble()

for(i in 1:rep){
  df <- tibble(known_b = 
                 rep(b,
                     each = n),
               n = n)
  x_15 <- rPLB(n = n, b = b[1])
  x_20 <- rPLB(n = n, b = b[2])
  x_25 <- rPLB(n = n, b = b[3])
  x_sample <- c(x_15, x_20, x_25)
  df$m <- x_sample
  df$rep <- i
  sim_out[[i]] <- df
}

sim_df <- bind_rows(sim_out)
(test_method <- sim_df %>%
    # mutate(date = as.Date(collectDate)) %>%
    group_by(rep, known_b) %>%
    #filter(siteID == "CARI") %>%
    # create list-column
    nest() %>% 
    mutate(method_compare = 
             map(data,
                 compare_slopes,
                 rsp_var = "m",
                 dw_range = dw_range)) %>%
    ungroup() %>%
    select(-data) %>%
    unnest(cols = method_compare))

test_method <- test_method %>%
  #filter(rep < 2) %>%
  mutate(in_CI = minCI < known_b & maxCI > known_b)

test_method %>%
  ggplot(aes(y = estimate, 
             ymin = minCI, 
             ymax = maxCI,
             x = rep,
             color = in_CI)) +
  geom_pointrange() +
  facet_wrap(~name)

test_method %>%
  group_by(name) %>%
  summarize(mean(in_CI))

test_method %>%
  ggplot(aes(x = known_b,
             y = estimate, 
             color = name)) +
  geom_point() +
  stat_smooth(method = "lm")

coef(lm(estimate ~ known_b, data = test_method %>%
          filter(name == "ELBn")))
coef(lm(estimate ~ known_b, data = test_method %>%
          filter(name == "MLE")))
coef(lm(estimate ~ known_b, data = test_method %>%
          filter(name == "NAS")))

test_method %>%
  group_by(known_b, name) %>%
  summarize(min = min(estimate),
            max = max(estimate)) %>%
  mutate(spread = max - min) %>%
  arrange(name)
