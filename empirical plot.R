# empirical data plot

# combine the AMD and NEON data into a single multipanel plot

library(ggplot2)
library(png)
library(grid)
# library(gridExtra)
library(ggpubr)
amd_a <- readPNG("figures/AMD_plot.png")
amd_b <- readPNG("figures/amd_relationship.png")
neon_a <- readPNG("figures/NEON_plot.png")
neon_b <- readPNG("figures/neon_relationship.png")
ggarrange(rasterGrob(amd_a),
             rasterGrob(amd_b),
             rasterGrob(neon_a),
             rasterGrob(neon_b),
             labels = c("A", "B", "C", "D"))

ggsave("figures/empirical_combined.png",
       height = 6,
       width = 4, 
       units = "in")
