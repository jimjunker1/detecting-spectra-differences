# sample size analysis

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")
library(tidybayes)

# read in datasets with variable sample size, n
# All other params are as in main analysis
n200 <- readRDS("data_sim/PLB_sim_n200_dat.rds")
n500 <- readRDS("data_sim/PLB_sim_n500_dat.rds")
n5000 <- readRDS("data_sim/PLB_sim_n5000_dat.rds")
n10000 <- readRDS("data_sim/PLB_sim_n10000_dat.rds")
n1000 <- readRDS("data_sim/PLB_sim_dat.rds")
n1000$n <- 1000


### ELBn struggles to estimate b when sampled from -2 and -2.5 PLB ###
# Filter out all reps which have an estimate == NA

# make vector of reps which have an NA value for the estimate
#rep_na <- n100$rep[which(is.na(n100$estimate))]

# filter out all the reps with NA values
# the !rep means (reps which do not match values in rep_na)
#PLB_sim_n100 <- n100 %>%
#  filter(!rep %in% rep_na)

# combine all data sets
dat <- bind_rows(n200, n500, n5000, n10000, n1000)

# labda ~n line graph old? ####
# # How does the lambda estimate vary by sample size?
# ggplot(dat, aes(x = n, y = estimate, color = name)) +
#   geom_point(alpha = 0.1) +
#   #stat_smooth(method = "loess") +
#   stat_summary(geom = "line",  fun = mean, size = 1.5) +
#   geom_hline(aes(yintercept = known_b), linetype = "dashed") +
#   facet_wrap(~known_b,
#              scales = "free") +
#   theme_bw() +
#   labs(title = "Lambda ~ n",
#        caption = "dashed line is known lambda value") +
#   NULL
# 
# ggsave("figures/lambda_n_5_sites.png",
#        width = 8,
#        height = 8)

# lambda ~ n density old? ####
# how does estimate of lambda vary with sample size?
# ggplot(dat,
#        aes(x = estimate,
#            y = ..scaled..,
#            fill = name)) +
#   geom_density(alpha = 0.5,
#                adjust = 2) +
#   facet_wrap(n ~ known_b,
#              ncol = 5,
#              labeller = label_both) +
#   geom_vline(aes(xintercept = known_b),
#              size = 1,
#              alpha = 0.75,
#              color = "black")+
#   scale_fill_viridis_d(option = "plasma") +
#   labs(title = "Sample size and known b",
#        x = "slope estimate") +
#   theme_bw() +
#   NULL

# lambda estimate halfeye ####
dat %>% 
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(
    aes(x = estimate, 
        y = Model,
        fill = Model)) +
  stat_halfeye(.width = c(0.66, 0.95)) +
  scale_fill_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  theme_bw() +
  geom_vline(aes(xintercept = known_b),
             linetype = "dashed") +
  labs(
    x = "Lambda estimate") +
  facet_wrap(n~known_b,
             ncol = 5,
             labeller = label_both,
             scales = "free_x") +
  theme(legend.position = "none") +
  NULL

ggsave(paste0("figures/", 
              substitute(n_vary),
              "_est_b.png"),
       width = 8,
       height = 8)

# regressions ####
# what are the estimated relationships across gradient with varying n?
# dat %>%
#   ggplot(aes(x = env_gradient,
#              y = estimate, 
#              group = rep,
#              color = rep)) +
#   stat_smooth(geom = "line",
#               method = "lm",
#               alpha = 0.15,
#               se = FALSE)+
#   geom_point() +
#   facet_wrap(n~name,
#              ncol = 3,
#              labeller = label_both)+
#   theme_bw()

dat %>% 
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(
    aes(x = env_gradient, 
        y = estimate, 
        color = name,
        group = rep)) +
  stat_smooth(method = "lm",
              se = FALSE, 
              geom = "line",
              size = 1.75,
              alpha = 0.15,
              color = "black") +
  geom_point() +
  facet_wrap(n~name,
             ncol = 3,
             labeller = label_both)+
  scale_color_manual(
    values = c("#FF914A", "#019AFF",  "#FF1984" )) +
  theme_bw()

ggsave(paste0("figures/", 
              substitute(n_vary),
              "_main.png"),
       width = 8,
       height = 8)

# relationship estimate halfeye ####
# distribution of relationship estimate (beta values) with varied n
relationship_estimate <- dat %>%
  group_by(rep, name, n) %>%
  nest() %>%
  mutate(lm_mod =
           map(data,
               ~lm(estimate ~ env_gradient, data = .x))) %>%
  mutate(tidied = map(lm_mod, broom::tidy)) %>%
  unnest(tidied) %>%
  filter(term == "env_gradient") %>%
  select(-data, -lm_mod, -statistic)

# old plot ####
# relationship_estimate %>%
#   ggplot(aes(y = ..scaled..,
#              x = estimate,
#              fill = name)) + 
#   geom_density(alpha = 0.5,
#                adjust = 2) +
#   geom_vline(xintercept = -0.5,
#              size = 1,
#              linetype = "dashed") +
#   theme_bw() +
#   scale_fill_viridis_d(option = "plasma") +
#   labs(x = "relationship estimate") +
#   facet_wrap(~n,
#              labeller = label_both,
#              ncol = 1) +
#   NULL
# new plot 8/30/22 ####
relationship_estimate %>%
  mutate(Model = factor(name,
                      levels = 
                        c("MLE",
                          "ELBn", 
                          "NAS"))) %>%
  ggplot(
    aes(x = estimate, 
        y = Model,
        fill = Model)) +
  stat_halfeye(.width = c(0.66, 0.95)) +
  scale_fill_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  theme_bw() +
  geom_vline(xintercept = -0.5,
                          size = 1,
                          linetype = "dashed")  +
  labs(
    x = "Relationship estimate") +
  facet_wrap(~n,
             scales = "free_x") +
  theme(legend.position = "none") +
  NULL

ggsave(filename = 
         paste0("figures/",
                substitute(n_vary),
                "_relationship_density.png"),
       width = 8,
       height = 8)
