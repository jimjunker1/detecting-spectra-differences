# script 1
# simulate data sets
# analyze using different methods
# Compare results


# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra)
source("custom_functions.R")


# Simulate ----------------------------------------------------------------


n <-  1000
b <- c(-1.5, -2, -2.5)
x_15 <- rPLB(n = n, b = b[1])
x_20 <- rPLB(n = n, b = b[2])
x_25 <- rPLB(n = n, b = b[3])
x_sample <- c(x_15, x_20, x_25)

df <- tibble(known_b = 
               rep(c(-1.5, -2, -2.5),
                   each = n),
             n = n,
             name = 
               rep(c("m_invert_1.5",
                     "m_invert_2.0",
                     "m_invert_2.5"),
                   each = n))
df$m <- x_sample

# just MLE estimate
# mle_lambda <- df %>%
#   group_by(known_b) %>%
#   nest() %>%
#   mutate(lambda = map(data,
#                       MLE_tidy,
#                       "m")) %>%
#   unnest(cols = lambda) %>%
#   select(-data) %>%
#   ungroup()

# compare slopes method from Pomeranz et al. 2022
dw_range = range(df$m, na.rm = TRUE)

(test_method <- df %>%
  # mutate(date = as.Date(collectDate)) %>%
  group_by(known_b) %>%
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

# Scale up simulation ####

n <-  1000
b <- c(-1.5, -2, -2.5)
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
