# sampling function

sites[sites$Sitio == 160,]

test <- sites |>
  st_drop_geometry() |>
  select(Sitio, Area, start, end) |>
  group_by(Sitio) |>
  arrange(start)

test$end <- test$end - 1

test$start <- rcarbon::BCADtoBP(test$start)
test$end <- rcarbon::BCADtoBP(test$end)
# remove duplicates
test <- test[!duplicated(test), ]

test$Area <- round(test$Area, 2)

test_lag <- sites |>
  group_by(Sitio) |>
  arrange(start) |> # if BP use desc()
  mutate(gap = start - lag(end)) |>
  select(Sitio, Fase, Area, start, end, gap) |>
  st_drop_geometry()


# sort(unique(test_lag$Sitio[test_lag$gap > 1]))

library(ggplot2)

ggplot(data = test_lag[test_lag$Sitio %in% sort(unique(test_lag$Sitio[test_lag$gap > 0])),]) +
  geom_line(aes(x = start, y = Area, group = Sitio),
            ) +
  geom_point(aes(x = start, y = Area), shape = 19, size = 2, color = 'red') +
  facet_wrap(~Sitio) +
  theme_minimal() +
  geom_vline(
    xintercept = sort(unique(test_lag$start)),
    linetype = "dashed",
    colour = "grey60"
  )



site_summary <- site_data %>%
  group_by(Sitio) %>%
  summarise(
    n_phases = first(n_phases),
    mean_duration = mean(duration),
    sd_duration = sd(duration),
    median_duration = median(duration),
    mean_growth_rate = mean(mean_growth_rate),
    sd_growth_rate = sd(mean_growth_rate),
    median_growth_rate = median(mean_growth_rate),
    q025_growth_rate = quantile(mean_growth_rate, 0.025),
    q975_growth_rate = quantile(mean_growth_rate, 0.975),
    .groups = 'drop'
  )

# Analyze interval-level results
interval_summary <- interval_data %>%
  group_by(Sitio, interval) %>%
  summarise(
    n_iterations = n(),
    mean_year_from = mean(year_from),
    mean_year_to = mean(year_to),
    mean_duration = mean(time_diff),
    mean_growth_rate = mean(growth_rate),
    sd_growth_rate = sd(growth_rate),
    median_growth_rate = median(growth_rate),
    q025_growth_rate = quantile(growth_rate, 0.025),
    q975_growth_rate = quantile(growth_rate, 0.975),
    .groups = 'drop'
  )

# View results
head(site_summary)
head(interval_summary)

# Example: Plot growth rate distribution for a specific site
library(ggplot2)

# Site-level growth rates
ggplot(site_data %>% filter(Sitio == 286), 
       aes(x = mean_growth_rate)) +
  geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7) +
  labs(title = "Distribution of Mean Growth Rates for Site 286",
       x = "Mean Growth Rate (ha/year)",
       y = "Count") +
  theme_minimal()

# Interval-level growth rates for multi-interval sites
ggplot(interval_data %>% filter(Sitio == 286), 
       aes(x = growth_rate, fill = factor(interval))) +
  geom_density(alpha = 0.5) +
  labs(title = "Growth Rate Distribution by Interval for Site 286",
       x = "Growth Rate (ha/year)",
       y = "Density",
       fill = "Interval") +
  theme_minimal()
