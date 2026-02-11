#### Functions used to compute longevity and hazard rates ----

fit_functions <- function(longevities, survival, set_func = 'all') {
  
  # Create survival object (all events are observed, so all get status = 1)
  surv_obj <- Surv(longevities, survival)
  
  func_names <- c('exponential', 'weibull', 'lognormal', 'gompertz', 'gamma', 'llogis', 'gengamma')
  
  results <- list()
  
  if (set_func == 'all') {
    for(func in func_names) {
      cat('Fitting', func, '...')
      results[[func]] <- flexsurvreg(surv_obj ~ 1, dist = func)
      cat(' ✓\n')
    }
  } else {
    
    # Validate input
    if (!all(set_func %in% func_names)) {
      stop(
        "Invalid distribution name(s).\n",
        "Allowed values are: 'all' or any of:\n  ",
        paste(func_names, collapse = ", ")
      )
    }
    
    # Fit selected
    for (func in set_func) {
      cat('Fitting', func, '...')
      results[[func]] <- flexsurvreg(surv_obj ~ 1, dist = func)
      cat(' ✓\n')
    }
  }
  
  return(results)
}


# ===== USAGE EXAMPLES =====

# Fit all distributions
all_fits <- fit_functions(times, set_func = 'all')

# Fit only specific ones
selected_fits <- fit_functions(times, set_func = c('exp'))

# Fit just two
two_fits <- fit_functions(times, set_func = c('exp', 'weib'))

# Access results
all_fits$weib$AIC
selected_fits$gengamma$coefficients

compare_functions <- function(longevities, survival, set_func = 'all') {
  
  function_list <- fit_functions(longevities, survival, set_func = 'all')
  
  comparison <- data.frame(
    Model = names(function_list),
    AIC = c(function_list[[1]]$AIC, function_list[[2]]$AIC, 
            function_list[[3]]$AIC, function_list[[4]]$AIC,
            function_list[[5]]$AIC, function_list[[6]]$AIC,
            function_list[[7]]$AIC),
    N_params = c(1, 2, 2, 2, 2, 2, 3)
  )
  
  # Sort by AIC and calculate Delta AIC and weights
  comparison <- comparison |>
    arrange(AIC) |>
    mutate(
      Delta_AIC = AIC - min(AIC),
      Weight = exp(-0.5 * Delta_AIC) / sum(exp(-0.5 * Delta_AIC)),
      Rank = row_number()
    )
  return(comparison)
}
