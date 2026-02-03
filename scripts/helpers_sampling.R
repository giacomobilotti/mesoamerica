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
    if (delay < delta) return(delay) # if <= is used, values that equal start/end could be sampled
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

## Convert BCAD to BP 
# I modified the two functions available in the rcarbon pacakge (https://cran.r-project.org/web/packages/rcarbon/)
# Since this is the only function I needed, I modified it according to my needs and stored locally.

bcad_to_bp <- function (x, bc_to_bp = TRUE) {
  index <- !is.na(x)
  if(bc_to_bp == TRUE) {
    if (any(x[index] == 0)) {
      stop("0 BC/AD is not a valid year.")
    }
    
    res <- matrix(c(x, rep(NA, length(x))), ncol = 2)
    res[index & x > 0, 2] <- abs(res[index & x > 0, 1] - 1950)
    res[index & x < 0, 2] <- abs(res[index & x < 0, 1] - 1949)
    return(res[, 2]) 
  } else if(bc_to_bp == FALSE) {
    
    if (any(x[index] < 0)) {
      stop("Post-bomb dates (<0 BP) are not currently supported.")
    }
    res <- matrix(c(x, rep(NA, length(x))), ncol = 2)
    res[index & x < 1950, 2] <- 1950 - res[index & x < 1950, 
                                           1]
    res[index & x >= 1950, 2] <- 1949 - res[index & x >= 1950, 
                                            1]
    return(res[, 2])
  }
}

#### Summarise MC sim results by site ----
summary_site <- function(sites) {
  
  tmp <- sites |>
    group_by(id, start, end) |>
    summarise(
      median     = median(growth_cagr),
      mean       = mean(growth_cagr),
      q025       = quantile(growth_cagr, 0.025, na.rm = TRUE),   
      q975       = quantile(growth_cagr, 0.975, na.rm = TRUE),   
      q25        = quantile(growth_cagr, 0.25, na.rm = TRUE),    
      q75        = quantile(growth_cagr, 0.75, na.rm = TRUE),    
      iqr        = IQR(growth_cagr, na.rm = TRUE), # Interquartile range
      # Linear growth
      med_lin    = median(growth_linear),
      q025_lin   = quantile(growth_linear, 0.025, na.rm = TRUE),
      q975_lin   = quantile(growth_linear, 0.975, na.rm = TRUE),
      # Year
      med_year   = median(bcad_to_bp(year, bc_to_bp = FALSE)),
      q025_year  = quantile(year, 0.025),
      q975_year  = quantile(year, 0.975),
      area       = mean(area) # the val is repeated the mean = the actual value 
    )
  return(tmp)
  
}

#### Plot ----
# a good function to create your own palettes with greyscale
# scales::alpha("steelblue", 0.2)
plot_site_cagr <- function(site, legend = TRUE, 
                           col_med = 'steelblue', col_mean = 'skyblue',
                           col_0 = 'grey10', col_95 = '#4682B433', col_50 = '#4682B44C') {
  
  required_cols <- c("median", "med_year", "q025", "q975", "q25", "q75")
  missing_cols <- setdiff(required_cols, colnames(site))
  
  if (length(missing_cols) > 0)
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  
  x_lim <- c(bcad_to_bp(max(site$start, na.rm = TRUE), bc_to_bp = FALSE),
             bcad_to_bp(min(site$end, na.rm = TRUE), bc_to_bp = FALSE))
  # Plot CAGR with SD
  plot(site$med_year, site$median,
       type = "b", pch = 16, col = NA,
       ylim = range(c(site$q025[-nrow(site)], site$q975[-nrow(site)]), na.rm = TRUE),
       xlim = x_lim,
       xlab = "Year (BC/AD)", ylab = "CAGR",
       main = paste("Growth Rate for Site", site$id[1]))
  
  # Add confidence interval polygon
  polygon(
    c(site$med_year[-nrow(site)], rev(site$med_year[-nrow(site)])),
    c(site$q025[-nrow(site)], rev(site$q975[-nrow(site)])),
    col = col_95, border = NA
  )
  polygon(
    c(site$med_year[-nrow(site)], rev(site$med_year[-nrow(site)])),
    c(site$q25[-nrow(site)], rev(site$q75[-nrow(site)])),
    col = col_50, border = NA
  )
  
  # Add mean line
  points(site$med_year, site$median, type = 'b', pch = 19,
         col = col_med, lwd = 2)
  lines(site$med_year, site$mean,
        col = col_mean, lwd = 2, lty = 2)
  
  # Add horizontal line at 0
  abline(h = 0, lty = 2, col = col_0)
  
  if(legend == TRUE) {
    legend("topleft", 
           legend = c("Median", "Mean", "50% CI", "95% CI"),
           col = c(col_med, col_mean, col_50, col_95),
           lty = c(1,2,1,1), lwd = c(2, 2, 10, 10), 
           pch = c(19, NA, NA, NA), bty = 'n') 
  }
}
# spaghetti plot function (valid for area and growth trajectories)
plot_spaghetti <- function(site, summary_data, iter_num = 50, seed = 1234, type = "growth_cagr",
                           median_col = "skyblue", x_lim = NULL, y_lim = NULL) {
  
  match.arg(type, choices = c("growth_cagr", "growth_linear", "area"))
  
  set.seed(seed)
  
  sample_iterations <- sample(unique(site$iteration), iter_num)
  site_sample <- site[site$iteration %in% sample_iterations, ]
  
  if(is.null(x_lim)) {
    # careful: END < START in years BP 
    x_lim <- c(bcad_to_bp(max(site$start, na.rm = TRUE), bc_to_bp = FALSE),
               bcad_to_bp(min(site$end, na.rm = TRUE), bc_to_bp = FALSE))
  }
  if(is.null(y_lim)) {
    y_lim <- c(
      min(site[[type]], na.rm = TRUE) - 0.1*min(site[[type]], na.rm = TRUE),
      max(site[[type]], na.rm = TRUE) + 0.1*max(site[[type]], na.rm = TRUE)
    )
  }
  
  cap_type <- strsplit(type, "_")[[1]][1]
  plot(NULL, xlim = x_lim, ylim = y_lim, xlab = "Year", 
       ylab = type, 
       main = paste0("Individual ", cap_type, " trajectories (", iter_num, " samples)")
  )
  
  for(iter in sample_iterations) {
    tmp <- site_sample[site_sample$iteration == iter, ]
    # for safety, but not necessary
    tmp <- tmp[order(tmp$year, decreasing = TRUE), ]
    lines(bcad_to_bp(tmp$year, bc_to_bp = FALSE), 
          tmp[[type]], col = rgb(0, 0, 0, 0.2))
    
  }
  
  if(type == "growth_cagr") {
    abline(h = 0, lty = 2, col = "gray50")
    # Add median on top
    lines(summary_data$med_year, summary_data$median, col = median_col, lwd = 3)
  } else if(type == "growth_linear") {
    abline(h = 0, lty = 2, col = "gray50")
    # Add median on top
    lines(summary_data$med_year, summary_data$med_lin, col = median_col, lwd = 3)
  } else if(type == "area") {
    lines(summary_data$med_year, summary_data$area, col = median_col, lwd = 3)
  }
  abline(v = bcad_to_bp(summary_data$start, bc_to_bp = FALSE), lty = 2, lwd =.75)
  
} 