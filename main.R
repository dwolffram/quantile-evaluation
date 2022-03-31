library(tidyverse)
source("plot_coverage.R")

# Load data
models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

# National level

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location == "US",
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

plot_coverage(df, B=100, type="confidence2", difference=FALSE)
ggsave("figures/national_coverage_confidence.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

plot_coverage(df, B=100, type="consistency", difference=FALSE)
ggsave("figures/national_coverage_consistency.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

# State level

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name != "US",
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

plot_coverage(df, B=100, type="confidence", difference=FALSE)
ggsave("figures/states_coverage.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

plot_coverage(df, B=100, type="confidence", difference=TRUE)
ggsave("figures/states_coverage_diff.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)


# Vermont

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name == "Vermont",
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

plot_coverage(df, B=100, type="confidence", difference=FALSE)
ggsave("figures/Vermont_coverage_confidence.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

plot_coverage(df, B=100, type="consistency", difference=FALSE)
ggsave("figures/Vermont_coverage_consistency.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)



### Multiple states at once

models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name %in% c("Nebraska", "Vermont", "Florida"),
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

coverage2 <- function(df) {
  df %>%
    group_by(model, location_name, quantile) %>%
    summarize(l = mean(truth < value), u = mean(truth <= value), .groups = "drop")
}

get_consistency_bands2 <- function(df) {
  df %>%
    group_by(model, quantile) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(get_consistency_interval(quantile, count, 0.9) / count,
           get_consistency_interval(quantile, count, 0.5) / count) %>%
    select(-count)
}

coverage_df <- coverage2(df)

bands <- get_consistency_bands2(df)

coverage_df <- coverage_df %>% 
  left_join(bands)

ggplot(coverage_df) +
  facet_grid(location_name ~ model) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
  geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
  geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
  geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
  scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1),
                     labels = function(x) ifelse(x == 0, "0", x)) +
  scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)) +
  xlab("Quantile level") +
  ylab("Coverage") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_line(size = 0.05)) +
  coord_fixed()

# ggsave("figures/coverage_3states_consistency.pdf", width=200, height=200, unit="mm", device = "pdf", dpi=300)