library(tidyverse)
library(tidybayes)

source("custom_functions.R")

# load data
est_lambda <- readRDS("data_sim/est_lambda_dat.rds") %>% 
  mutate(target = known_b,
         target_name = "Lambda") #figure 1 sims

angles = readRDS("data_sim/angles.rds") %>% 
  mutate(target = known_relationship,
         minCI = estimate - 1.96*std.error,
         maxCI = estimate + 1.96*std.error,
         target_name = "Regression Slope") #figure 3 sims

# rel_data = readRDS("data_sim/rel_data.rds") # figure 4 sims

rel_data_summary = readRDS("data_sim/rel_data_summary.rds") %>% 
  mutate(target = known_relationship,
         minCI = estimate - 1.96*std.error,
         maxCI = estimate + 1.96*std.error,
         target_name = "Regression Slope") # figure 5 sims

# estimate bias
lambda_bias = est_lambda %>% 
  mutate(diff = estimate - known_b) %>% 
  group_by(name, known_b) %>% 
  mutate(reps = max(rep)) %>%
  summarize(bias = sum(diff>0)/max(reps),
            diff_mean = mean(diff),
            diff_sd = sd(diff)) %>% 
  mutate(figure = "Figure 1",
         target = "Lambda")

angle_bias = angles %>% 
  mutate(diff = estimate - known_relationship) %>% 
  group_by(id, name, known_relationship) %>% 
  mutate(reps = max(rep)) %>%
  summarize(bias = sum(diff>0)/max(reps),
            diff_mean = mean(diff),
            diff_sd = sd(diff)) %>% 
  mutate(figure = "Figure 3",
         target = "Regression Slope")


slope_bias = rel_data_summary %>% 
  mutate(diff = estimate - known_relationship) %>% 
  group_by(name, known_relationship, term) %>% 
  mutate(reps = max(rep)) %>%
  summarize(bias = sum(diff>0)/max(reps),
            diff_mean = mean(diff),
            diff_sd = sd(diff)) %>% 
  mutate(figure = "Figure 5",
         target = "Regression Slope")



#plot and summarize
  
all_dat = bind_rows(est_lambda, rel_data_summary, angles) %>% 
  filter(!is.na(minCI)) %>% 
  mutate(target_zero = target - target,
         estimate_z = estimate - target,
         upper_z = maxCI - target,
         lower_z = minCI - target) %>%
  arrange(-estimate_z) %>% 
  mutate(rank = 1:nrow(.))


all_dat  %>%
  group_by(name,target_name) %>% 
  arrange(estimate_z) %>% 
  mutate(rank = rank(estimate_z)) %>% 
  filter(rep <= 5) %>% 
  ggplot(aes(y = rank, x = estimate_z, xmin = lower_z, xmax = upper_z, 
             color = name)) + 
  geom_pointrange(size = 0.1) + 
  geom_vline(xintercept = 0) + 
  facet_wrap(target_name~name, scales = "free_y")


bias_table = all_dat %>% 
  mutate(conf_width = maxCI - minCI,
         diff = estimate - target,
         abs_bias = abs(diff)) %>% 
  group_by(target_name, name) %>%
  add_count %>% 
  group_by(target_name, name, n) %>%
  summarize(median_ci_range = median(conf_width),
            median_abs_bias = median(abs_bias),
            sd_abs_bias = sd(abs_bias)) %>% 
  arrange(target_name, median_ci_range)

write_csv(bias_table, file = "tables/bias_table.csv")


all_dat %>% 
  mutate(conf_width = maxCI - minCI,
         diff = estimate - target,
         abs_bias = abs(diff)) %>%
  group_by(target_name) %>%
  mutate(mean = mean(diff),
         sd = sd(diff),
         diff_s = (diff - mean)/sd) %>% 
  ggplot(aes(x = diff_s)) + 
  geom_histogram() + 
  facet_wrap(target_name ~ name) +
  # scale_x_log10() +
  NULL
  

