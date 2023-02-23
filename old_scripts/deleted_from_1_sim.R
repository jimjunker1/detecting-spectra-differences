# old stuff removed from 1_simulations.R
# PLB no relationship
# for reproducibility
set.seed(2806)

PLB_static_b <- sim_result(n = 1000,
                           b = rep(-2, length(env_gradient)),
                           env_gradient = env_gradient,
                           rep = rep, 
                           m_lower = m_lower,
                           m_upper = m_upper,
                           distribution = "PLB")
saveRDS(PLB_static_b, file = paste0("data_sim/",
                                    substitute(PLB_static_b),
                                    "_dat.rds"))
plot_sim(PLB_static_b)
rm(PLB_static_b)


# PLB "small" relationship
# lambda varies from -1.9 to -2.1
# for reproducibility
set.seed(2806)

PLB_sim_small <- sim_result(n = 1000,
                            b = seq(from = -1.9, to = -2.1, length.out = 5),
                            env_gradient = env_gradient,
                            rep = rep, 
                            m_lower = m_lower,
                            m_upper = m_upper,
                            distribution = "PLB")

saveRDS(PLB_sim_small, file = paste0("data_sim/",
                                     substitute(PLB_sim_small),
                                     "_dat.rds"))

plot_sim(PLB_sim_small)
rm(PLB_sim_small)

# sim n200

# re run ------------------------------------------------------------------


# Re run this ####
# threw an error in group 517
# tomorrow, attempt to use purrr::possibly() as a wrapper
# should hopefully throw an NA instead of stopping the run...
# https://purrr.tidyverse.org/reference/safely.html 

# seems to throw lots or errors or have too few bins
# changing to n = 200 on 6/13

## manuscript --> sample size not the main objective, but anecdotally seems that n= 100 is too few observations. 



# addedd purrr::possibly() to custom functions
# double check this simulation and make sure it looks right
# also double check the number of iterations it thre out (only one window popped up while I was watching...)
# make sure plot_sim() function still works with NA's... may return a list???
# may need to modify code in subsequent calculations
# i.e. does the number of reps actually = 1000???

# PLB_sim_n100 %>% drop_na(estimate) %>%
#   plot_sim()
# 
# PLB_sim_n100 %>% drop_na(estimate) %>%
#   group_by(name) %>%
#   count()

#rep_na <- PLB_sim_n100$rep[which(is.na(PLB_sim_n100$estimate))]

# PLB_sim_n100 %>%
#   filter(!rep %in% rep_na, !is.na(name)) %>%
#   select(-method_compare) %>%
#   drop_na() %>%
#   plot_sim()

# sim_data <- PLB_sim_n100
# 
# relationship_estimate <- sim_data %>%
#   filter(!rep %in% rep_na, !is.na(name)) %>%
#   select(-method_compare) %>%
#   group_by(rep, name, known_relationship) %>%
#   nest() %>%
#   mutate(lm_mod =
#            map(data,
#                ~lm(estimate ~ env_gradient, data = .x))) %>%
#   mutate(tidied = map(lm_mod, broom::tidy)) %>%
#   unnest(tidied) %>%
#   filter(term == "env_gradient") %>%
#   select(-data, -lm_mod, -statistic)
# 
# slope_distribution <- relationship_estimate %>%
#   group_by(name) %>%
#   summarize(mu_slope = mean(estimate, na.rm = TRUE),
#             sd_slope = sd(estimate),
#             p25 = quantile(estimate, probs = 0.25),
#             p50 = quantile(estimate, probs = 0.5),
#             p975 = quantile(estimate, probs = 0.975))


# Shallow lambda ####
# lambda parameter varies between -1.1 and -1.5 across a hypothetical gradient from -1 to 1. i.e. relationship of beta = 0.2
# setting main analysis to 5 sites
lambda_shallow = seq(from = -1.1, to = -1.5, length.out = 5)

# for reproducibility
set.seed(2806)

PLB_shallow_lambda <- sim_result(n = 1000,
                                 b = lambda_shallow,
                                 env_gradient = env_gradient,
                                 rep = rep, 
                                 m_lower = m_lower,
                                 m_upper = m_upper,
                                 distribution = "PLB")

saveRDS(PLB_shallow_lambda, file = paste0("data_sim/",
                                          substitute(PLB_shallow_lambda),
                                          "_dat.rds"))
plot_sim(PLB_shallow_lambda)