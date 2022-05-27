# Other distributions

# size spectra are "expected" to follow a bounded power law distribution (Anderson, Edwards et al.), but it is not know what distribution they will follow in reality. 

# this script samples from alternate distributions in order to test the robustness of the results sampling from a bounded power law. 

# distributions
# pareto, truncated pareto, weibull


# Preliminaries -----------------------------------------------------------

library(tidyverse)
library(sizeSpectra)
#library(EnvStats)
library(VGAM)
source("custom_functions.R")


# truncated pareto ####
# b = shape
# tpareto shape parameter has to be positive
# High shape values (i.e. shape = 2) make it unlikely to sample large sizes
# rescaled for comparison
# Might drop tpareto sample. May confuse the issue too much
rep = 1000
tp_sim <- sim_result(b = c(0.1, 0.6, 1.1),
                     env_gradient = c(-1, 0, 1),
                     rep = rep, 
                     m_lower = 0.0026,
                     m_upper = 1.2e4,
                     distribution = "tpareto")

plot_sim(tp_sim)
rm(tp_sim)