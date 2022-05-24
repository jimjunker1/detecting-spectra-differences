# custom functions
# these are largely copied from Pomeranz et al. 2018 Freshwater Biology and Pomeranz et al. 2022 Global Change Biology
# copied 5/24/22 --> go through these and make sure doing what I think I remember them doing. 

# MLE in tidy format
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
                  p = PL.bMLE, x = x, n = length(x), 
                  xmin = xmin, xmax = xmax,
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
  AS_lm <- lm(log_count~log_mids_center,
              data = NAS)
  AS_out <- est_ci(AS_lm)
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
  ELB_lm <- lm(log_count~log_mids_center,
             data = ELB)
  # ELB_coef <- coef(ELB_lm)
  # ELB_conf <- confint(ELB_lm)
  ELB_out <- est_ci(ELB_lm)
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
