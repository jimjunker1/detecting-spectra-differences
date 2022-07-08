# Main simulation


## this takes about 4 hours to run...
# I'm sure the time could be improved with the use of parallel computing
# But I have a 4 month old and let the sims run over night instead of going down an optimizing/troubleshooting rabbit hole. 
# probably less time efficient in the long run but {shruggy guy emoji}


## added "possibly()" command to custom functions
# load PLB_sim, save as something else and rerun simulation
# Does old (pre-possibly()) data match new data with possibly()?
# Need to re-run all sims with modified functions??

## Possibly didn't seem to work, chnaged n=100 to n=200 and it seems to be good. Need to go through and delete comments / clean it all up. 

# for reproducibility
set.seed(2806)

# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")

# Simulation and summary --------------------------------------------------

# PLB ####
# lambda parameter varies between -1.5 and -2.5 across a hypothetical gradient from -1 to 1. i.e. relationship of beta = 0.5
# setting main analysis to 5 sites
lambda = seq(from = -1.5, to = -2.5, length.out = 5)
env_gradient = seq(from = -1, to = 1, length.out = 5)

lambda
env_gradient

rep = 1000
# body sizes in NEON data set varies from 0.0026 to 1.2e4
# setting sample range based on that minimum, but reducing max by order of magnitude
m_lower = 0.0026
m_upper = 1.2e3

debugonce(sim_result)
PLB_sim <- sim_result(n = 1000,
                      b = lambda,
                      env_gradient = env_gradient,
                      rep = rep, 
                      m_lower = m_lower,
                      m_upper = m_upper,
                      distribution = "PLB")

saveRDS(PLB_sim, file = paste0("data_sim/",
                               substitute(PLB_sim),
                               "_dat.rds"))

PLB_sim = readRDS(file = paste0("data_sim/",
                                substitute(PLB_sim),
                                "_dat.rds"))
plot_sim(PLB_sim)
#rm(PLB_sim)

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

# sensitivity, supplementary information ----------------------------------

# different scale of env_gradient (predictor, x)
# for reproducibility
set.seed(2806)

PLB_large_x <- sim_result(n = 1000,
                          b = lambda,
                          env_gradient = seq(from = -100,
                                             to = 100,
                                             length.out = 5),
                          rep = rep, 
                          m_lower = m_lower,
                          m_upper = m_upper,
                          distribution = "PLB")

saveRDS(PLB_large_x, file = paste0("data_sim/",
                               substitute(PLB_large_x),
                               "_dat.rds"))

plot_sim(PLB_large_x)
rm(PLB_large_x)

# change range of m sampled in rPLB
# for reproducibility
set.seed(2806)

PLB_small_m <- sim_result(n = 1000,
                          b = lambda,
                          env_gradient = env_gradient,
                          rep = rep, 
                          m_lower = 1,
                          m_upper = 100,
                          distribution = "PLB")

saveRDS(PLB_small_m, file = paste0("data_sim/",
                               substitute(PLB_small_m),
                               "_dat.rds"))

plot_sim(PLB_small_m)
rm(PLB_small_m)

# change number of sites across "gradient"
# 3 sites
# for reproducibility
set.seed(2806)

PLB_3_sites <- sim_result(
  b = c(-1.5, -2, -2.5),
  env_gradient = c(-1, 0, 1),
  rep = rep, 
  m_lower = m_lower,
  m_upper = m_upper,
  distribution = "PLB")

saveRDS(PLB_3_sites, file = paste0("data_sim/",
                               substitute(PLB_3_sites),
                               "_dat.rds"))

plot_sim(PLB_3_sites)
rm(PLB_3_sites)

# 10 sites
# for reproducibility
set.seed(2806)

PLB_10_sites <- sim_result(
  b = round(seq(from = -1.5, to = -2.5, length.out = 10), 2),
  env_gradient = round(seq(from = -1, to = 1, length.out = 10), 2),
  rep = rep, 
  m_lower = m_lower,
  m_upper = m_upper,
  distribution = "PLB")

saveRDS(PLB_10_sites, file = paste0("data_sim/",
                                   substitute(PLB_10_sites),
                                   "_dat.rds"))

plot_sim(PLB_10_sites)
rm(PLB_10_sites)



# sample size of body sizes -----------------------------------------------
# varying the sample size of simulated data
# 10,000, 5,000, 500, 100

#rep = 1000
n = 10000 
# for reproducibility
set.seed(2806)

PLB_sim_n10000 <- sim_result(n = n,
                             b = lambda,
                             env_gradient = env_gradient,
                             rep = rep, 
                             m_lower = m_lower,
                             m_upper = m_upper,
                             distribution = "PLB")
PLB_sim_n10000$n <- n

plot_sim(PLB_sim_n10000)

saveRDS(PLB_sim_n10000, file = paste0("data_sim/",
                                    substitute(PLB_sim_n10000),
                                    "_dat.rds"))

#rm(PLB_sim_n10000)

n = 5000 
# for reproducibility
set.seed(2806)

PLB_sim_n5000 <- sim_result(n = n,
                             b = lambda,
                             env_gradient = env_gradient,
                             rep = rep, 
                             m_lower = m_lower,
                             m_upper = m_upper,
                             distribution = "PLB")

plot_sim(PLB_sim_n5000)

PLB_sim_n5000$n <- n

saveRDS(PLB_sim_n5000, file = paste0("data_sim/",
                                      substitute(PLB_sim_n5000),
                                      "_dat.rds"))

#rm(PLB_sim_n5000)

n = 500 
# for reproducibility
set.seed(2806)

PLB_sim_n500 <- sim_result(n = n,
                            b = lambda,
                            env_gradient = env_gradient,
                            rep = rep, 
                            m_lower = m_lower,
                            m_upper = m_upper,
                            distribution = "PLB")

plot_sim(PLB_sim_n500)

PLB_sim_n500$n <- n

saveRDS(PLB_sim_n500, file = paste0("data_sim/",
                                     substitute(PLB_sim_n500),
                                     "_dat.rds"))

#rm(PLB_sim_n500)



# re run ------------------------------------------------------------------


# Re run this ####
# threw an error in group 517
# tomorrow, attempt to use purrr::possibly() as a wrapper
# should hopefully throw an NA instead of stopping the run...
# https://purrr.tidyverse.org/reference/safely.html 

# seems to throw lots or errors or have too few bins
# changing to n = 200 on 6/13

## manuscript --> sample size not the main objective, but anecdotally seems that n= 100 is too few observations. 

n = 200 
# for reproducibility
set.seed(2806)
PLB_sim_n200 <- sim_result(n = n,
                           b = lambda,
                           env_gradient = env_gradient,
                           rep = rep, 
                           m_lower = m_lower,
                           m_upper = m_upper,
                           distribution = "PLB")

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

plot_sim(PLB_sim_n200)

PLB_sim_n200$n <- n

saveRDS(PLB_sim_n200, file = paste0("data_sim/",
                                    substitute(PLB_sim_n200),
                                    "_dat.rds"))
rm(PLB_sim_n200)

# end re run --------------------------------------------------------------



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
