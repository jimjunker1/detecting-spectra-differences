# main plot version 2

# if we like this replace in custom functions script and re run main simulation. 

library(tidyverse)

sim_data <- readRDS("data_sim/PLB_sim_dat.rds")
sim_data %>%
  mutate(Model = factor(name,
                        levels = 
                          c("MLE",
                            "ELBn", 
                            "NAS"))) %>%
  ggplot(aes(x = env_gradient,
             y = estimate, 
             group = rep,
             color = Model)) +
  scale_color_manual(
    values = c("#019AFF", "#FF914A", "#FF1984" )) +
  stat_smooth(geom = "line",
              method = "lm",
              alpha = 0.15,
              se = FALSE, 
              color = "black") +
  geom_point(alpha = 0.5) +
  facet_wrap(~Model)+
  labs(title = sim_data$distribution,
       subtitle = 
         paste0("m_range =(",
                sim_data$m_lower,
                ",",
                sim_data$m_upper,
                ") \nb = (",
                min(sim_data$known_b),
                ",",
                max(sim_data$known_b),
                ")")) +
  theme_bw() +
  theme(legend.position="none")
ggsave("figures/PLB_main_v2.png", width = 2100, height = 2100, units = "px")
