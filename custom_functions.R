# custom functions
# these are largely copied from Pomeranz et al. 2018 Freshwater Biology and Pomeranz et al. 2022 Global Change Biology
# copied 5/24/22 --> go through these and make sure doing what I think I remember them doing.

# These functions require the tidyverse (Wickham et al.) and sizeSpectra (Edwards et al. 207, 2020)

# MLE in tidy format
# this is a wrapper function I wrote to be able to use the MLE functions from size Spectra, but used on data in the tidy format. 
# it was designed to specifically work with list-columns of the data you want to analyze. 
# I make no guarantees on it's functionality 

MLE_tidy <- function(df, rsp_var){
  # define variables
  x <- df[[rsp_var]]
  xmin = min(x)
  xmax = max(x)
  log.x = log(x)
  sum.log.x = sum(log.x)
  
  # initial starting point for parameter estimate
  PL.bMLE = 1/(log(min(x)) - sum.log.x/length(x)) - 1
  
  # non-linear minimization  
  PLB.minLL = nlm(negLL.PLB, 
                  p = PL.bMLE,
                  x = x,
                  n = length(x), 
                  xmin = xmin,
                  xmax = xmax,
                  sumlogx = sum.log.x)
  
  # estimate for b
  PLB.bMLE = PLB.minLL$estimate
  # minimum estimate of b
  PLB.minNegLL = PLB.minLL$minimum
  
  ## 95% CI calculation
  bvec = seq(PLB.bMLE - 0.5, PLB.bMLE + 0.5, 1e-05)
  PLB.LLvals = vector(length = length(bvec))
  for (i in 1:length(bvec)) {
    PLB.LLvals[i] = negLL.PLB(bvec[i],
                              x = x,
                              n = length(x), 
                              xmin = xmin,
                              xmax = xmax,
                              sumlogx = sum.log.x)
  }
  critVal = PLB.minNegLL + qchisq(0.95, 1)/2
  bIn95 = bvec[PLB.LLvals < critVal]
  # confidence interval
  PLB.MLE.bConf = c(min(bIn95), max(bIn95))
  if (PLB.MLE.bConf[1] == min(bvec) | 
      PLB.MLE.bConf[2] == max(bvec)) {
    dev.new()
    plot(bvec, PLB.LLvals)
    abline(h = critVal, col = "red")
    stop("Need to make bvec larger - see R window")
  }
  # return b estimate and min/max 95% CI
  return(data.frame(b = PLB.bMLE,
                    minCI = min(bIn95),
                    maxCI = max(bIn95)))
}
# binning method
bin_and_center <- function(data, var, breaks, ...){
  # data is a data frame
  # var is a string, and is the name of a column in data which you want to bin
  # breaks controls the number of bins as defined in hist(). Can be a vector giving breakpoints of bins [i.e. Log10 breaks = 10^seq(min, max)], a single number designating number of bins, etc. See ?hist for details. If not supplied by the user a default of log2 width bins are calculated 
  
  # are breaks supplied by the user?
  if(exists("breaks") == FALSE){
    
    # calculate log2 width bins which are inclusive of the range of data supplied
    breaks = 2^seq(floor(range(log2(data[[var]]))[1]),
                   ceiling(range(log2(data[[var]]))[2]))
    message("breaks not supplied, using log2 width bins")
  }
  
  # bin values using hist()
  binned_hist = hist(data[[var]], 
                     breaks = breaks,
                     include.lowest = TRUE, plot = FALSE)
  # calculate "left" and "right" edge of bins
  breaks_orig = binned_hist$breaks[1:(length(breaks)-1)]
  breaks_offset = binned_hist$breaks[2:length(breaks)]
  # total bin width = right edge - left edge
  bin_width = breaks_offset - breaks_orig
  count = binned_hist$counts
  log_mids = log10(binned_hist$mids)
  biomass = count * 10**log_mids
  nbiomass = log10(biomass / bin_width)
  dataout = data.frame(
    count = count,
    log_count = log10(count),
    # normalize counts =count / width (White et al 1997)
    log_count_corrected = log10(count / bin_width),
    # original midpoint of bin log10 transformed
    log_mids = log_mids,
    bin_width = bin_width,
    biomass = biomass,
    nbiomass = nbiomass)
  # remove bins with 0 counts
  # -Inf comes from log10(count / break_width) above
  dataout = dataout[dataout$log_count_corrected !=-Inf,]
  # # recenter data at x=0
  # mid_row = ceiling(nrow(dataout)/2)
  # # subtract value of mid row from all mids
  # dataout$log_mids_center = 
  #   dataout[,"log_mids"] - dataout[mid_row,"log_mids"]
  dataout
}
# function to extract slope estimate and CI's
est_ci <- function(lm_obj){
  out <- c(coef(lm_obj)[2],
    confint(lm_obj)[2],
    confint(lm_obj)[4])
  names(out) <- c("estimate", "minCI", "maxCI")
  out
}

# function to compare slopes####
compare_slopes <- function(data, dw_range, rsp_var, ...){
  # define breaks widths and mid bins for log2 and 6 equal log bins
  # log2
  breaks2 <- 2^(floor(log2(min(dw_range))):
                  ceiling(log2(max(dw_range))) )
  mid_bin_2 <- log10(breaks2[floor(length(breaks2)/2)]) 
  # 6 equal log bins
  breaks_log_6 <- exp(seq(floor(log(min(dw_range))),
                          ceiling(log(max(dw_range))),
                          length.out = 7))
  mid_bin_log_6 <- log10(breaks_log_6[floor(length(breaks_log_6)/2)])
  
  # Normalized Abundance Spectra
  NAS <- bin_and_center(data, var = rsp_var, breaks = breaks2)
  NAS$log_mids_center <- NAS$log_mids - mid_bin_2
  NAS_lm <- lm(log_count_corrected~log_mids_center,
               data = NAS)
  NAS_out <- est_ci(NAS_lm)
  # NAS_coefs <- coef(NAS_lm)
  # NAS_conf <- confint(NAS_lm)
  # non-normalized abundance
  # AS_lm <- lm(log_count~log_mids_center,
  #             data = NAS)
  # AS_out <- est_ci(AS_lm)
  # AS_coefs<- coef(AS_lm)
  # AS_conf <- confint(AS_lm)
  
  # equal logarithmic binning method
  
  # breaks_log_6 <- exp(seq(floor(log(min(dw_range))),
  #                         ceiling(log(max(dw_range))),
  #                         length.out = 7))
  # mid_bin_log_6 <- log10(breaks_log_6[floor(length(breaks_log_6)/2)])
  # ELB6
   ELB <- bin_and_center(data, var = rsp_var, breaks = breaks_log_6)
  ELB$log_mids_center <- ELB$log_mids - mid_bin_log_6
  # ELB_lm <- lm(log_count~log_mids_center,
  #            data = ELB)
  # # ELB_coef <- coef(ELB_lm)
  # # ELB_conf <- confint(ELB_lm)
  # ELB_out <- est_ci(ELB_lm)
  ELBn_lm <- lm(log_count_corrected~log_mids_center,
              data = ELB)
  # ELBn_coef <- coef(ELBn_lm)
  # ELBn_conf <- confint(ELBn_lm)
  ELBn_out <- est_ci(ELBn_lm)
  
  # MLE method from Edwards et al. 2018
  mle_b <- MLE_tidy(data, rsp_var)
  names(mle_b) <- c("estimate", "minCI", "maxCI")
  
  slopes <- as.data.frame(rbind(NAS_out,
                   #AS_out,
                   #ELB_out,
                   ELBn_out,
                   mle_b))
  slopes$name <- c("NAS", #"AS", 
                   #"ELB",
                   "ELBn", "MLE")
  # slopes = data.frame(NAS = NAS_out,
  #                     # NAS_25 = NAS_conf[2],
  #                     # NAS_975 = NAS_conf[4],
  #                     AS = AS_out,
  #                     # AS_25 = AS_conf[2],
  #                     # AS_975 = AS_conf[4],
  #                     ELB6 = ELB_out,
  #                     # ELB_25 = ELB_conf[2],
  #                     # ELB_95 = ELB_conf[4],
  #                     ELBn = ELBn_out,
  #                     # ELBn_25 = ELBn_conf[2],
  #                     # ELBn_975 = ELBn_conf[4],
  #                     mle = mle_b
  #                     # mle_25 = mle_b$minCI,
  #                     # mle_975 = mle_b$maxCI
  #                     )
  slopes
}


# simulation funciton -----------------------------------------------------

# function to automate simulation runs

# Written 5/25/22

sim_result <- function(n = 1000,
                       b,
                       env_gradient,
                       rep,
                       m_lower,
                       m_upper,
                       distribution){
  stopifnot(length(b) == length(env_gradient))
  known_relationship <- -(max(b) - min(b)) / 
    (max(env_gradient) - min(env_gradient)) 
  # simulate values from distribution
  sim_out <- list()
  for(i in 1:rep){
    df <- tibble(
      known_b = rep(b, each = n),
      known_relationship = known_relationship,
      env_gradient = rep(env_gradient, each = n),
      n = n,
      distribution = distribution,
      m_lower = m_lower,
      m_upper = m_upper,)
    if(distribution == "PLB"){
      sample_list <- list()
      for(n_b in 1:length(b)){
        sample_list[[n_b]] <- rPLB(n = n,
                                  xmin = m_lower,
                                  xmax = m_upper,
                                  b = b[n_b])
      }
      m_sample <- unlist(sample_list)
}
    if(distribution == "tpareto"){
      sample_list <- list()
      for(n_b in 1:length(b)){
        sample_list[[n_b]] <- rtruncpareto(n = n,
                                   lower = m_lower,
                                   upper = m_upper,
                                   shape = b[n_b])
      }
      m_sample <- unlist(sample_list)
    }
    df$m <- m_sample
    df$rep <- i
    df$dist = distribution
    sim_out[[i]] <- df
  }
  
  sim_df <- bind_rows(sim_out)
  
  dw_range = range(sim_df$m, na.rm = TRUE)
  
  out <- sim_df %>%
    group_by(rep,
             known_b,
             known_relationship,
             env_gradient,
             distribution,
             m_lower,
             m_upper) %>%
    nest() %>% 
    mutate(method_compare = 
             map(data, 
                 # added the possibly() function on 6/12/2022
                 # small sample size (n=100) was having a hard time with MLE estimates for steep (b = -2.5) distributions
                 # added this so that it "skips" problematic iterations instead of stopping. 
                 # removed 6/13 - work this out later
                 compare_slopes,
                 rsp_var = "m",
                 dw_range = dw_range)) %>%
    ungroup() %>%
    select(-data) %>%
    unnest(cols = method_compare)
  out
}


plot_sim <- function(sim_data,
                     display = FALSE,
                     adjust = 1){
  # plot each regression "rep"
  main_plot <- ggplot(sim_data,
         aes(x = env_gradient,
             y = estimate, 
             group = rep,
             color = rep)) +
    stat_smooth(geom = "line",
                method = "lm",
                alpha = 0.15,
                se = FALSE)+
    geom_point() +
    facet_wrap(~name)+
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
    theme_bw()
  
  ggsave(plot = main_plot,
         filename = 
           paste0("figures/",
                  substitute(sim_data),
                  "_main.png"))
  
if(sim_data$distribution[1] == "PLB"){
  #### this is the old plot, keeping for now ####
  # estimate_density <- sim_data %>%
  #   ggplot(aes(y = ..scaled..,
  #              x = estimate,
  #              fill = name)) +
  #   geom_density(alpha = 0.5,
  #                adjust = adjust) +
  #   theme_bw() +
  #   facet_wrap(~known_b) +
  #   geom_vline(aes(xintercept = known_b),
  #              size = 1,
  #              linetype = "dashed",
  #              color = "black")+
  #   scale_fill_viridis_d(option = "plasma") +
  #   labs(title = sim_data$distribution,
  #        x = "slope estimate") +
  #   NULL
  #### New 8/30/22 "halfeye" plot
  estimate_density <- sim_data %>%
    mutate(Model = factor(name,
                          levels = 
                            c("MLE",
                              "ELBn", 
                              "NAS"))) %>%
    ggplot(
      aes(x = estimate, 
          y = Model,
          fill = Model)) +
    stat_halfeye(.width = c(0.66, 0.95)) +
    scale_fill_manual(
      values = c("#019AFF", "#FF914A", "#FF1984" )) +
    theme_bw() +
    geom_vline(aes(xintercept = known_b),
               linetype = "dashed") +
    labs(
      x = "Lambda estimate") +
    facet_wrap(~known_b,
               scales = "free_x")
  }
  
  if(sim_data$distribution[1] == "tpareto"){
    estimate_density <- sim_data %>%
      mutate(trans_b = known_b * -1 -1) %>%
      ggplot(aes(y = ..scaled..,
                 x = estimate,
                 fill = name)) +
      geom_density(alpha = 0.5,
                   adjust = adjust) +
      theme_bw() +
      facet_wrap(~trans_b) +
      geom_vline(aes(xintercept = (known_b*-1)-1),
                 size = 1,
                 linetype = "dashed",
                 color = "black")+
      scale_fill_viridis_d(option = "plasma") +
      labs(title = sim_data$distribution,
           subtitle = "(Shape parameter * -1) -1",
           x = "slope estimate") +
      NULL}
  
  ggsave(plot = estimate_density,
         filename = 
           paste0("figures/",
                  substitute(sim_data),
                  "_est_b_density.png"))
  
  # distribution of slope estimates?
  relationship_estimate <- sim_data %>%
    group_by(rep, name, known_relationship) %>%
    nest() %>%
    mutate(lm_mod =
             map(data,
                 ~lm(estimate ~ env_gradient, data = .x))) %>%
    mutate(tidied = map(lm_mod, broom::tidy)) %>%
    unnest(tidied) %>%
    filter(term == "env_gradient") %>%
    select(-data, -lm_mod, -statistic)
  
  slope_distribution <- relationship_estimate %>%
    group_by(name) %>%
    summarize(mu_slope = mean(estimate, na.rm = TRUE),
              sd_slope = sd(estimate, na.rm = TRUE),
              p25 = quantile(estimate, probs = 0.25, na.rm = TRUE),
              p50 = quantile(estimate, probs = 0.5, na.rm = TRUE),
              p975 = quantile(estimate, probs = 0.975, na.rm = TRUE))
  
  write_csv(slope_distribution,
            file = paste0("results/",
                          substitute(sim_data),
                          ".csv"))
  # plot density of relationship estimate
  # old relationship distribution plot ####  
  # relationship_density <- relationship_estimate %>%
    # ggplot(aes(y = ..scaled..,
    #            x = estimate, 
    #            fill = name)) + 
    # geom_density(alpha = 0.5,
    #              adjust = adjust) +
    # geom_vline(aes(xintercept = known_relationship),
    #            linetype = "dashed",
    #            size = 1)+
    # theme_bw() +
    # scale_fill_viridis_d(option = "plasma") +
    # labs(x = "relationship estimate") +
    # NULL
  # new 8/30/22 relationship halfeye plot ####
  relationship_density <- relationship_estimate %>%
    mutate(Model = factor(name,
                          levels = 
                            c("MLE",
                              "ELBn", 
                              "NAS"))) %>%
    ggplot(aes(x = estimate, 
               y = Model,
               fill = Model))+
    stat_halfeye(.width = c(0.66, 0.95)) +
    scale_fill_manual(
      values = c("#019AFF",
                 "#FF914A",
                 "#FF1984" )) +
    theme_bw() +
    geom_vline(
      aes(xintercept = known_relationship),
               linetype = "dashed") +
    labs(x = "Relationship estimate") +
    NULL
  
  ggsave(plot = relationship_density,
         filename = 
           paste0("figures/",
                  substitute(sim_data),
                  "_relationship_density.png"))
  
  if(display == TRUE){
    return(
      list(
        slope_distribution,
        main_plot,
        estimate_density,
        relationship_density))
  }
}

