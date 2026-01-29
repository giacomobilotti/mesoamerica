##### helper functions for the sampling

exp_sampling <- function(start, end, perc) {
  
  if (start - end < 0) stop("start - end is negative! Dates must be in years BP.")
  if (perc <= 0) stop("perc must be > 0")
  
  # phase length in years, always positive
  delta <- start - end
  
  # calculate lambda
  L <- 1 / (delta / perc) 
  
  for (i in 1:1000) {
    delay <- floor(rexp(1, rate = L))
    if (delay <= delta) return(delay)
  }
  stop("Failed to sample truncated exponential after 1000 attempts.\nIf 'perc' is very low, increase it; otherwise change seed / try again.")
}


norm_peak <- function(start, end, peak_sd_perc) {
  
  if (start <= end) stop("Start must be older than end (years BP).")
  if(peak_sd_perc < 0) stop("The peak_sd_perc you chose is invalid! It must be > 0.")
  if(peak_sd_perc > 1) warning("The peak_sd_perc you chose is > 1. This will make your sd bigger than the phase interval.")
  
  phase <- start - end 
  
  mean_peak <- start - phase / 2
  
  if (peak_sd_perc == 0) {
    warning("sd will be 0. All values equal the mean (no variation).")
    return(floor(mean_peak))
  }
  
  std <- phase*peak_sd_perc
  
  for(i in 1:1000){
    peak <- floor(rnorm(n = 1, mean = mean_peak, sd = std))
    if(peak <= start && peak >= end) return(peak)
  } 
  stop("Failed to sample truncated normal after 1000 attempts.\nIf 'sd' is large, decrease it; otherwise change seed / try again.")
  
}

check_sampled_dates <- function(df) {
  
  if (any(df$year <= 0)) warning("Non-positive years detected")
  
  by(df, df$id, function(x) {
    
    if (!isTRUE(all(diff(x$year) < 0))) {
      stop("Years are not strictly decreasing for site ", x$id[1])
    }
    
    if (any(x$area < 0)) {
      stop("Negative area detected for site ", x$id[1])
    }
  })
  
  invisible(TRUE)
}
