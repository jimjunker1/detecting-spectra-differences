# elegant analysis



# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra) #bounded power law and MLE functions
library(tidybayes) # for nicer looking plots
source("custom_functions.R")

# Simulation and summary --------------------------------------------------

# 1 lambda estimates ####
# lambda parameter varies between -2.5 and -0.5 
# env_gradient is used to identify sites, no regression estimate is performed here

# for reproducibility
set.seed(2806)

# 9 sites for facet_wrap aesthetics
lambda = seq(from = -0.5, to = -2.5, length.out = 9)
env_gradient = seq(from = -1, to = 1, length.out = 9)

lambda
env_gradient

rep = 1000
# body sizes in NEON data set varies from 0.0026 to 1.2e4
# setting sample range based on that minimum, but reducing max by order of magnitude
m_lower = 0.0026
m_upper = 1.2e3

est_lambda <- sim_result(n = 1000,
                      b = lambda,
                      env_gradient = env_gradient,
                      rep = rep, 
                      m_lower = m_lower,
                      m_upper = m_upper,
                      distribution = "PLB")

plot_sim(est_lambda)

# match colors, jitter points, increase alpha, try GAM/loess, 
  
(known_est_b_line <- est_lambda %>%
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(
       aes(x = known_b, 
           y = estimate, 
           color = name)) +
  geom_point(
    position = position_jitter(w = 0.025, h = 0),
    alpha = 0.25
  ) +
  stat_smooth(method = "lm",
              se = FALSE, 
              geom = "line",
              size = 1.75) +
  geom_abline(linetype = "dashed",
              size = 1) +
  scale_color_manual(
    values = c("#FF914A", "#019AFF",  "#FF1984" )) +
  theme_bw())

ggsave(plot = known_est_b_line,
       filename = "figures/known_est_b_line.png")

saveRDS(est_lambda, file = paste0("data_sim/",
                               substitute(est_lambda),
                               "_dat.rds"))

# 2 lambda windows ####

# 2A ####
# steep lambda [-1.5, -2.5]
# for reproducibility
set.seed(2806)

lambda = seq(from = -1.5, to = -2.5, length.out = 5)
env_gradient = seq(from = -1, to = 1, length.out = 5)

lambda
env_gradient

steep_lambda <- sim_result(n = 1000,
                         b = lambda,
                         env_gradient = env_gradient,
                         rep = rep, 
                         m_lower = m_lower,
                         m_upper = m_upper,
                         distribution = "PLB")

plot_sim(steep_lambda)

saveRDS(steep_lambda, file = paste0("data_sim/",
                                  substitute(steep_lambda),
                                  "_dat.rds"))

# 2b ####
# medium lambda [-1, -2]
# for reproducibility
set.seed(2806)

lambda = seq(from = -1, to = -2, length.out = 5)
env_gradient = seq(from = -1, to = 1, length.out = 5)

lambda
env_gradient

med_lambda <- sim_result(n = 1000,
                           b = lambda,
                           env_gradient = env_gradient,
                           rep = rep, 
                           m_lower = m_lower,
                           m_upper = m_upper,
                           distribution = "PLB")

plot_sim(med_lambda)

saveRDS(med_lambda, file = paste0("data_sim/",
                                    substitute(med_lambda),
                                    "_dat.rds"))

# 2c ####
# shallow lambda [-0.5, -1.5]
# for reproducibility
set.seed(2806)

lambda = seq(from = -0.5, to = -1.5, length.out = 5)

lambda
env_gradient

shallow_lambda <- sim_result(n = 1000,
                           b = lambda,
                           env_gradient = env_gradient,
                           rep = rep, 
                           m_lower = m_lower,
                           m_upper = m_upper,
                           distribution = "PLB")

plot_sim(shallow_lambda)

saveRDS(shallow_lambda, file = paste0("data_sim/",
                                  substitute(shallow_lambda),
                                  "_dat.rds"))

# 3 changing known relationship ####

# 1 env_gradient [-1, 1]
# 3 sets of lambdas [-1.5], [-1.25, -1.75], [-1, -2]

# 3A ####
# lambda = [-1.5]

set.seed(2806)

lambda = seq(from = -1.5, to = -1.5, length.out = 5)

lambda
env_gradient

relationship_0 <- sim_result(n = 1000,
                             b = lambda,
                             env_gradient = env_gradient,
                             rep = rep, 
                             m_lower = m_lower,
                             m_upper = m_upper,
                             distribution = "PLB")

plot_sim(relationship_0 )

saveRDS(relationship_0,
        file = paste0("data_sim/",
                      substitute(relationship_0),
                      "_dat.rds"))


# 3B ####
# lambda = [-1.25, -1.75]

set.seed(2806)

lambda = seq(from = -1.25, to = -1.75, length.out = 5)

lambda
env_gradient

relationship_025 <- sim_result(n = 1000,
                             b = lambda,
                             env_gradient = env_gradient,
                             rep = rep, 
                             m_lower = m_lower,
                             m_upper = m_upper,
                             distribution = "PLB")

plot_sim(relationship_025 )

saveRDS(relationship_025,
        file = paste0("data_sim/",
                      substitute(relationship_025),
                      "_dat.rds"))

# 3C ####
# lambda = [-1, -2]

set.seed(2806)

lambda = seq(from = -1, to = -2, length.out = 5)

lambda
env_gradient

relationship_05 <- sim_result(n = 1000,
                               b = lambda,
                               env_gradient = env_gradient,
                               rep = rep, 
                               m_lower = m_lower,
                               m_upper = m_upper,
                               distribution = "PLB")

plot_sim(relationship_05)

saveRDS(relationship_05,
        file = paste0("data_sim/",
                      substitute(relationship_05),
                      "_dat.rds"))
