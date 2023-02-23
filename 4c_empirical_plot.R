# empirical data plot

# combine the AMD and NEON data into a single multipanel plot

library(ggplot2)
library(ggpubr)
amd_a <- readRDS("figures/AMD_plot.rds")
amd_b <- readRDS("figures/amd_relationship.rds")
neon_a <- readRDS("figures/NEON_plot.rds")
neon_b <- readRDS("figures/neon_relationship.rds")

ggarrange(amd_a,
             amd_b,
             neon_a,
             neon_b,
             labels = c("A", "B", "C", "D"))

ggsave("figures/empirical_combined.png",
       height = 8,
       width = 8, 
       units = "in")
