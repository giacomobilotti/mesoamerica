#### Longevity  ----
## Use framework from https://www.pnas.org/doi/10.1073/pnas.2218834120
## It seems we can assume sites age like humans. This can be tested and the hazard/saturation can be calculated

# # load MC sim results (or run the script)
# all_sites_mc <- readRDS(file.path(targetdir, 'mc_unifnorm.rds'))

# calculate persistence 
mc_pers <- mc_persistence(all_sites_mc)

## add size information
mc_pers$size[mc_pers$id %in% large_sites] <- 'Large'
mc_pers$size[mc_pers$id %in% ex_large] <- 'Failed Large'
mc_pers$size[is.na(mc_pers$size)] <- 'Small'

## add survival information 
mc_pers$survival <- ifelse(
  mc_pers$norm_pers == 1, 0, 1
  # 0 = alive, 1 = dead
)

## Plot persistence (or longevity)
# gg_long <- 
  ggplot(mc_pers[mc_pers$norm_pers < 1,], aes(x = persistence)) +
  geom_histogram(aes(y = after_stat(density)), 
                 bins = 40, fill = "grey70", colour = "white") +
  theme_minimal() +
  labs(
    x = "Persistence (years)",
    y = "Count",
    title = "Site persistence"
  )
ggsave(filename = 'site_persistence_MC.tiff', gg_long, dpi = 1200, width = 7, height = 7, units = 'in')

# # # Facet by size category
# ggplot(mc_pers, aes(x = persistence, fill = size)) +
#   geom_histogram(
#     bins = 40,
#     position = "identity",
#     alpha = 0.5,
#     colour = "white"
#   ) +
#   theme_minimal() +
#   labs(
#     x = "Persistence (years)",
#     y = "Density",
#     title = "Site persistence by size class"
#   ) +
#   geom_density(alpha = 0.3, adjust = 1.2, aes(color = size)) +
#   facet_wrap(~size)
 
# ggplot(mc_pers, aes(x = persistence, colour = size, fill = size)) +
#   geom_density(alpha = 0.3, adjust = 1.2) +
#   theme_minimal()

# visually, it seems that the best fitting function could be a saturating one (try the density plots)
# this aligns with the results of pnas.2218834120
# However, we can also perform the same tests they did to find the best fitting function

# Load required packages
library(tidyverse)
library(survival)
library(flexsurv)  # This package handles multiple parametric survival models
library(MASS)      # For fitdistr

# check for the best fitting function for each iteration

# in 997 out of 1000 cases the best distribution is the generalised gamma dist.
# This makes sense as many distribution used for survival analysis (e.g., Gamma, exp, Weibull)
# are special cases of this one
# best_f <- lapply(
#   X = 1:max(mc_pers$iteration),
#   FUN = function(i) {
#     times <- mc_pers$persistence[mc_pers$iteration == i]
#     results <- compare_functions(times) 
#     return(results)
#   }
# )
# 
# bests <- c()
# second <- c()
# for(i in seq_along(best_f)) {
#   bests <- c(bests, best_f[[i]]$Model[1])
#   second <- c(second, best_f[[i]]$Model[2])
# }

all_fun <- lapply(
  X = 1:max(mc_pers$iteration),
  FUN = function(i) {
    times <- mc_pers$persistence[mc_pers$iteration == i]
    results <- fit_functions(times)
    return(results)
  }
)


i <- 10
times <- mc_pers$persistence[mc_pers$iteration == i]

hist(times, breaks = 50, probability = TRUE)

pars <- all_fun[[i]]$gengamma$res

curve(
  dgengamma(x,
            mu = pars["mu","est"],
            sigma = pars["sigma","est"],
            Q = pars["Q","est"]),
  add = TRUE,
  col = "red",
  lwd = 2
)

hist(
  rgengamma(
    n = 1000,
    mu = pars["mu","est"],
    sigma = pars["sigma","est"],
    Q = pars["Q","est"])
  )
)


mus <-  c()
sigmas <- c()
Qs <- c()
shapes <- c()
rates <- c()

for(i in 1:1000) {
  coefs <- all_fun[[i]]$gengamma$coefficients
  # shapes <- c(shapes, coefs[1])
  # rates <- c(rates, coefs[2])
  mus <- c(mus, coefs[1])
  sigmas <- c(sigmas, coefs[2])
  Qs <- c(Qs, coefs[3])
}

all_fun[[1]]$gamma$res

coefs <- data.frame(
  mu = mus,
  sigma = sigmas,
  Q = Qs #,
  # shape = shapes,
  # rate = rates
)

best_c <- colMeans(coefs)

for(i in 1:1000){
  hist(rgengamma(
  n = 100,
  shape = coefs$mu[i],
  sigma = coefs$sigma[i],
  Q = coefs$Q[i]
  ))
  }


hist(rgamma(
  n = 1000,
  shape = exp(best_c[1]),
  rate  = exp(best_c[2])
))

all_fun[[1]]$gengamma$coefficients
?rgengamma()
mean(mc_pers$persistence[mc_pers$norm_pers < 1], na.rm = TRUE)

# Get the actual shape parameter (not the log-transformed one)
weibull_shape <- exp(fit_weib$coefficients["shape"])
weibull_scale <- exp(fit_weib$coefficients["scale"])

cat("Weibull shape:", weibull_shape, "\n")
cat("Weibull scale:", weibull_scale, "\n")

# Interpretation:
# shape < 1: decreasing hazard (early failures)
# shape = 1: constant hazard (exponential)
# shape > 1: increasing hazard (wear out)


#####


t_grid <- seq(min(times), max(times), length.out = 500)

weib_par <- fit_weib$res[, "est"]
shape <- weib_par["shape"]
scale <- weib_par["scale"]

weib_density <- dweibull(t_grid, shape = shape, scale = scale)

ggplot(data.frame(times = times), aes(x = times)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 40,
                 fill = "grey80",
                 colour = "white") +
  theme_minimal() +
  labs(
    x = "Persistence (years)",
    y = "Density",
    title = "Persistence distribution with fitted survival models"
  ) + geom_line(
    data = data.frame(t = t_grid, d = weib_density),
    aes(x = t, y = d),
    colour = "red",
    linewidth = 1
  )

library(survival)

km <- survfit(surv_obj ~ 1)

plot(
  km,
  xlab = "Persistence (years)",
  ylab = "Survival probability",
  conf.int = FALSE,
  col = "black",
  lwd = 2
)

S_weib <- pweibull(t_grid, shape = shape, scale = scale, lower.tail = FALSE)

lines(t_grid, S_weib, col = "red", lwd = 2)

ln_par <- fit_lnorm$res[, "est"]

S_lnorm <- plnorm(
  t_grid,
  meanlog = ln_par["meanlog"],
  sdlog   = ln_par["sdlog"],
  lower.tail = FALSE
)

lines(t_grid, S_lnorm, col = "blue", lwd = 2)

h_weib <- (shape / scale) * (t_grid / scale)^(shape - 1)

plot(
  t_grid, h_weib,
  type = "l",
  lwd = 2,
  col = "red",
  xlab = "Persistence (years)",
  ylab = "Hazard"
)

#####
t_grid <- seq(min(times), max(times), length.out = 500)

ll_par <- fit_llogis$res[, "est"]
shape <- ll_par["shape"]
scale <- ll_par["scale"]

h_llogis <- (shape / scale) *
  (t_grid / scale)^(shape - 1) /
  (1 + (t_grid / scale)^shape)

w_par <- fit_weib$res[, "est"]
h_weib <- (w_par["shape"] / w_par["scale"]) *
  (t_grid / w_par["scale"])^(w_par["shape"] - 1)

plot(
  t_grid, h_llogis,
  type = "l",
  lwd = 2,
  col = "red",
  ylim = c(0, max(h_llogis, h_weib)),
  xlab = "Persistence (years)",
  ylab = "Hazard"
)

lines(t_grid, h_weib, col = "blue", lwd = 2, lty = 2)

legend(
  "topright",
  legend = c("Log-logistic (saturating)", "Weibull (decreasing)"),
  col = c("red", "blue"),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)

S_llogis <- 1 / (1 + (t_grid / scale)^shape)


plot(
  km,
  conf.int = FALSE,
  lwd = 2,
  xlab = "Persistence (years)",
  ylab = "Survival probability"
)

lines(t_grid, S_llogis, col = "red", lwd = 2)

gg_par <- fit_gengamma$res[, "est"]
mu    <- gg_par["mu"]
sigma <- gg_par["sigma"]
Q     <- gg_par["Q"]

library(flexsurv)

h_gg <- hgengamma(t_grid, mu = mu, sigma = sigma, Q = Q)

plot(
  t_grid, h_gg,
  type = "l",
  lwd = 2,
  col = "black",
  xlab = "Persistence (years)",
  ylab = "Hazard",
  main = "Implied hazard (generalised gamma)"
)
lines(t_grid, h_llogis, col = "red", lwd = 2, lty = 2)
legend(
  "topright",
  legend = c("Generalised gamma", "Log-logistic"),
  col = c("black", "red"),
  lwd = 2,
  lty = c(1, 2),
  bty = "n"
)



mc_pers2 <- mc_pers |>
  filter(
    is.finite(persistence),
    persistence > 0
  )

mc_pers2 <- mc_pers2 |>
  mutate(persistence_sc = persistence / 100)

library(flexsurv)
library(survival)

fit_gengamma_safe <- function(df) {
  
  surv_obj <- Surv(df$persistence_sc, rep(1, nrow(df)))
  
  flexsurvreg(
    surv_obj ~ 1,
    dist = "gengamma",
    inits = c(
      mu    = mean(log(df$persistence_sc)),
      sigma = sd(log(df$persistence_sc)),
      Q     = 0.5
    ),
    control = list(maxit = 2000)
  )
}

fits_gengamma <- by(
  mc_pers2,
  mc_pers2$size,
  fit_gengamma_safe
)

ggplot(mc_pers2, aes(x = persistence, fill = size)) +
  geom_histogram(
    aes(y = after_stat(density)),
    bins = 40,
    alpha = 0.4,
    position = "identity"
  ) +
  theme_minimal()

tgrid <- seq(
  0,
  quantile(mc_pers2$persistence_sc, 0.99),
  length.out = 400
)

dens_df <- do.call(rbind, lapply(names(fits_gengamma), function(cl) {
  
  fit <- fits_gengamma[[cl]]
  pars <- fit$res
  
  mu    <- pars["mu","est"]
  sigma <- pars["sigma","est"]
  Q     <- pars["Q","est"]
  
  dens <- dgengamma(
    x = tgrid,
    mu = mu,
    sigma = sigma,
    Q = Q
  )
  
  data.frame(
    persistence = tgrid * 100,
    density = dens,
    size = cl
  )
}))


ggplot(mc_pers2, aes(x = persistence, fill = size)) +
  geom_histogram(
    aes(y = 150*after_stat(density)),
    bins = 40,
    alpha = 0.35,
    position = "identity"
  ) +
  geom_line(
    data = dens_df,
    aes(x = persistence, y = density, colour = size),
    linewidth = 1.2
  ) +
  theme_minimal()

# sites still occupied at the terminus
still_alive <- mc_pers2 %>%
  filter(end == 431)   # or whatever represents "still active"

simulate_abandonment <- function(n_sims = 1000, t0 = 0, size_class = "Large") {
  
  # Get the fitted parameters for the class
  fit <- fits_gengamma[[size_class]]
  pars <- fit$res
  mu    <- pars["mu","est"]
  sigma <- pars["sigma","est"]
  Q     <- pars["Q","est"]
  
  # Simulate from the fitted gengamma
  sims <- rgengamma(n_sims, mu = mu, sigma = sigma, Q = Q)
  
  # Add current age of each site
  sims <- sims + t0  # t0 = years already survived
  return(sims)
}

set.seed(123)

simulated_times <- still_alive %>%
  rowwise() %>%
  mutate(
    sim_abandon = list(
      simulate_abandonment(
        n_sims = 1,
        t0 = persistence,   # current age
        size_class = size
      )
    )
  ) %>%
  unnest(cols = c(sim_abandon))

n_sims <- 100

simulated_times_mc <- still_alive %>%
  rowwise() %>%
  mutate(
    sim_abandon = list(
      simulate_abandonment(
        n_sims = n_sims,
        t0 = persistence,
        size_class = size
      )
    )
  ) %>%
  unnest(cols = c(sim_abandon))

ggplot(simulated_times_mc, aes(x = sim_abandon, fill = size)) +
  geom_histogram(bins = 40, alpha = 0.4, position = "identity") +
  theme_minimal() +
  labs(
    x = "Simulated total persistence (years)",
    y = "Count",
    title = "Predicted abandonment times for still-occupied sites"
  )

