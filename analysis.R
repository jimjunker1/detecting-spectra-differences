# Cleaned up simulation

# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
source("custom_functions.R")


# Simulation and summary --------------------------------------------------

# PLB ####
rep = 1000
m_lower = 0.0026
m_upper = 1.2e4

start_reg = proc.time()
PLB_sim <- sim_result(b = c(-1.5, -2, -2.5),
                      env_gradient = c(-1, 0, 1),
                      rep = rep, 
                      m_lower = m_lower,
                      m_upper = m_upper,
                      distribution = "PLB")
plot_sim(PLB_sim)
rm(PLB_sim)

# PLB no relationship
PLB_static_b <- sim_result(b = c(-2, -2, -2),
                      env_gradient = c(-1, 0, 1),
                      rep = rep, 
                      m_lower = m_lower,
                      m_upper = m_upper,
                      distribution = "PLB")
plot_sim(PLB_static_b)
rm(PLB_static_b)

proc.time() - start_reg

# PLB "small" relationship
PLB_sim_small <- sim_result(b = c(-1.9, -2, -2.1),
                      env_gradient = c(-1, 0, 1),
                      rep = rep, 
                      m_lower = m_lower,
                      m_upper = m_upper,
                      distribution = "PLB")
plot_sim(PLB_sim_small)
rm(PLB_sim_small)

# sensitivity, supplementary information ----------------------------------

# different scale of env_gradient (predictor, x)
PLB_large_x <- sim_result(b = c(-1.5, -2, -2.5),
                          env_gradient = c(-100, 0, 100),
                          rep = rep, 
                          m_lower = m_lower,
                          m_upper = m_upper,
                          distribution = "PLB")

plot_sim(PLB_large_x)
rm(PLB_large_x)

# change range of m sampled in rPLB
PLB_small_m <- sim_result(b = c(-1.5, -2, -2.5),
                          env_gradient = c(-1, 0, 1),
                          rep = rep, 
                          m_lower = 1,
                          m_upper = 100,
                          distribution = "PLB")

plot_sim(PLB_small_m)
rm(PLB_small_m)

# add more sites across "gradient"
PLB_10_sites <- sim_result(
  b = round(seq(from = -1.5,
                to = -2.5,
                length.out = 10), 2),
  env_gradient = round(seq(from = -1,
                     to = 1,
                     length.out = 10), 2),
  rep = rep, 
  m_lower = m_lower,
  m_upper = m_upper,
  distribution = "PLB")

plot_sim(PLB_10_sites)
rm(PLB_10_sites)


