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
  

# proportion of 95% CI with target ####
lambda_ci_prop <- est_lambda %>%
  select(known_b, rep, estimate, minCI, maxCI, name) %>%
  mutate(in_ci = known_b > minCI & known_b < maxCI) %>%
  na.omit() %>%
  group_by(name) %>%
  summarize(count = n(),
            proportion = sum(in_ci, na.rm = TRUE) / count) %>%
  mutate(target_name = "Lambda")

beta_ci_prop <- angles %>%
  select(known_relationship, rep, estimate, minCI, maxCI, name, id) %>%
  mutate(in_ci = known_relationship > minCI & known_relationship< maxCI) %>%
  na.omit() %>%
  group_by(name) %>%
  summarize(count = n(),
            proportion = sum(in_ci, na.rm = TRUE) / count) %>%
  mutate(target_name = "lambda scenarios")

rel_beta_ci_prop <- rel_data_summary %>%
  select(known_relationship, rep, estimate, minCI, maxCI, name) %>%
  mutate(in_ci = known_relationship > minCI & known_relationship< maxCI) %>%
  na.omit() %>%
  group_by(name) %>%
  summarize(count = n(),
            proportion = sum(in_ci, na.rm = TRUE) / count) %>%
  mutate(target_name = "Varying beta")

beta_props <- bind_rows(beta_ci_prop, rel_beta_ci_prop) %>%
  group_by(name) %>%
  summarise(count = sum(count), 
            proportion = mean(proportion)) %>%
  mutate(target_name = "Regression Slope")



ci_prop <- bind_rows(lambda_ci_prop, beta_props)

table <- left_join(bias_table, ci_prop)
write_csv(table, file = "tables/bias_table.csv")

# plot of 95% CI's ####
# estimated lambda
est_lambda %>%
  filter(rep < 500) %>%
  select(known_b, rep, estimate, minCI, maxCI, name) %>%
  group_by(name, known_b) %>%
  mutate(in_ci = known_b > minCI & known_b < maxCI,
         rank = rank(minCI)) %>%
  na.omit() %>%
  ggplot(aes(x = estimate, 
             xmin = minCI, 
             xmax = maxCI,
             color = in_ci,
             y = rank)) +
  geom_linerange() +
  facet_wrap(known_b~name,
             scales = "free") +
  theme_bw()

ggsave("figures/lamba_ci_true.png",
       width = 8,
       height = 8)

rel_data_summary %>%
  filter(rep < 500) %>%
  select(known_relationship, rep, estimate, minCI, maxCI, name) %>%
  group_by(name, known_relationship) %>%
  mutate(in_ci = known_relationship > minCI & known_relationship< maxCI,
         rank = rank(minCI)) %>%
  na.omit() %>%
  ggplot(aes(x = estimate, 
             xmin = minCI, 
             xmax = maxCI,
             color = in_ci,
             y = rank)) +
  geom_linerange() +
  facet_wrap(known_relationship~name,
             scales = "free") +
  theme_bw()

ggsave("figures/rel_data_ci_true.png",
       width = 8,
       height = 8)

angles %>%
  filter(rep < 500) %>%
  select(known_relationship, rep, estimate, minCI, maxCI, name, id) %>%
  group_by(name, known_relationship) %>%
  mutate(in_ci = known_relationship > minCI &
           known_relationship < maxCI,
         rank = rank(minCI)) %>%
  na.omit() %>%
  ggplot(aes(xmin = minCI, 
             xmax = maxCI,
             color = in_ci,
             y = rank)) +
  geom_linerange() +
  facet_wrap(id~name,
             scales = "free") +
  theme_bw()

ggsave("figures/angles_ci_true.png",
       width = 8,
       height = 8)
