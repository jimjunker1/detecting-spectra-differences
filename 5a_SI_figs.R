# figures and tables for SI

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")
library(tidybayes)

# body values sample size ####

# read in datasets with variable sample size, n
# All other params are as in main analysis
n200 <- readRDS("data_sim/PLB_sim_n200_dat.rds")
n500 <- readRDS("data_sim/PLB_sim_n500_dat.rds")
n5000 <- readRDS("data_sim/PLB_sim_n5000_dat.rds")
n10000 <- readRDS("data_sim/PLB_sim_n10000_dat.rds")
# n1000 are the results from the main analysis
# the actual number of body sizes sampled is n = 999
n1000 <- readRDS("data_sim/PLB_sim_dat.rds")
n1000$n <- 999


# combine all sample size data sets data sets ####
dat <- bind_rows(n200, n500, n5000, n10000, n1000)


## absolute deviation based on body sample size ####
# modify data to calulate absolute difference in estimate and known lamba
abs_dev_n_vary <- dat %>% 
  mutate(diff_abs = abs(estimate - known_b))

# summary table ####
# write table of estimates based on sample size
abs_dev_n_vary %>% 
  mutate(body_n = n) %>%
  group_by(name, body_n) %>%
  summarize(diff_abs_mean = mean(diff_abs, na.rm = TRUE),
             diff_abs_sd = sd(diff_abs, na.rm = TRUE)) %>%
  write_csv("tables/abs_dev_n_vary.csv")

# summary plot ####
# plot of absolute deviation
abs_dev_n_vary %>%
  ggplot(aes(x = n, y = diff_abs, color = name)) +
  geom_point(alpha = 0.05,
             position = position_dodge(200)) +
  stat_summary(fun=mean, 
               geom="line", 
               size=1.5) +
  #geom_line() +
  theme_bw() +
  scale_color_manual(
    values = c("#FF1984",
               "#FF914A",
               "#019AFF"),
    breaks = c("L2n",
               "ELBn", 
               "MLE")) +
  labs(x = "number of sampled body sizes",
       y = "absolute deviation of estimate")

ggsave("figures/abs_dev_by_body_n.png",
       width = 8,
       height = 8)


# range of M and envs_gradient ####

# Changing range of m, range of hypothetical gradient
main_result <- n1000
PLB_large_x  <- readRDS("data_sim/PLB_large_x_dat.rds")
PLB_small_m <- readRDS("data_sim/PLB_small_m_dat.rds")

main_result$id <- "main"
PLB_large_x$id <- "large_x"
PLB_small_m$id <- "small_m"


# summary table SI 2 ####
# bind data by rows
si_2_dat <- bind_rows(main_result,
                      PLB_large_x,
                      PLB_small_m) 
# calc and write summary table
si_2_dat %>%
  mutate(diff_abs = abs(estimate - known_b)) %>%
  group_by(id, name) %>%
  summarize(diff_abs_mean = mean(diff_abs, na.rm = TRUE),
            diff_abs_sd = sd(diff_abs, na.rm = TRUE)) %>%
  arrange(name, id) %>%
  select(name, id, diff_abs_mean, diff_abs_sd) %>%
  write_csv("tables/si_table2.csv")

# number of sites ####
PLB_3_sites <- readRDS("data_sim/PLB_3_sites_dat.rds")
PLB_10_sites <- readRDS("data_sim/PLB_10_sites_dat.rds")
PLB_3_sites$id <- "sites_03"
PLB_10_sites$id <- "sites_10"
main_result$id <- "sites_05"

# summary table si 3 ###
# bind data by rows
si_3_dat <- bind_rows(main_result,
                      PLB_3_sites,
                      PLB_10_sites) 
# calc and write summary table
si_3_dat %>%
  mutate(diff_abs = abs(estimate - known_b)) %>%
  group_by(id, name) %>%
  summarize(diff_abs_mean = mean(diff_abs, na.rm = TRUE),
            diff_abs_sd = sd(diff_abs, na.rm = TRUE)) %>%
  arrange(name, id) %>%
  select(name, id, diff_abs_mean, diff_abs_sd) %>%
  write_csv("tables/si_table3.csv")
