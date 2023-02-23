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
library(tidybayes) # for nicer looking plots
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
plot_sim(PLB_sim)
#rm(PLB_sim)





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

#PLB_small_m <- readRDS("data_sim/PLB_small_m_dat.rds")

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

plot_sim(PLB_sim_n200)

PLB_sim_n200$n <- n

saveRDS(PLB_sim_n200, file = paste0("data_sim/",
                                    substitute(PLB_sim_n200),
                                    "_dat.rds"))
rm(PLB_sim_n200)

