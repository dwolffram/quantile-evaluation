library(tidyverse)
library(patchwork)
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

p1 <- plot_coverage(df, B=100, type="consistency", difference=FALSE) +
  facet_wrap("Consistency") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank())
p2 <- plot_coverage(df, B=100, type="confidence2", difference=FALSE)+
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.background.y = "grey",
        strip.text)

p1 / p2

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

### Multiple types

models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

df1 <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name == "US",
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

df2 <- df1

df1$type <- "consistency"
df2$type <- "confidence"

df <- bind_rows(df1, df2)

coverage2 <- function(df) {
  df %>%
    group_by(model, type, quantile) %>%
    summarize(l = mean(truth < value), u = mean(truth <= value), .groups = "drop")
}

get_consistency_bands2 <- function(df) {
  df %>%
    group_by(model, type, quantile) %>%
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



library(tidyverse)
Sys.setlocale("LC_ALL", "C")

coverage <- function(df) {
  df %>%
    group_by(model, quantile) %>%
    summarize(l = mean(truth < value), u = mean(truth <= value), .groups = "drop")
}

# resample groups
# in our case dates - with multiple rows for different locations, models, etc.
sample_n_groups = function(grouped_df, n, replace = FALSE) {
  groups_resampled <- grouped_df %>% 
    group_keys() %>% 
    slice_sample(n = n, replace = replace)
  
  grouped_df %>% 
    right_join(groups_resampled, by = group_vars(grouped_df))
}

get_confidence_bands <- function(df, date_column, B = 100) {
  coverage_df <- data.frame()
  
  for (i in 1:B) {
    df_resampled <- df %>%
      group_by(!!date_column) %>%
      sample_n_groups(n = n_distinct(group_keys(.)), replace = TRUE)
    
    coverage_sample <- coverage(df_resampled)
    
    coverage_df <- bind_rows(coverage_df, coverage_sample)
  }
  
  # compute CIs from bootstrapped coverage
  coverage_df %>%
    group_by(model, quantile) %>%
    summarize(l_5 = quantile(l, 0.05),
              l_95 = quantile(l, 0.95),
              l_25 = quantile(l, 0.25),
              l_75 = quantile(l, 0.75),
              u_5 = quantile(u, 0.05),
              u_95 = quantile(u, 0.95),
              u_25 = quantile(u, 0.25),
              u_75 = quantile(u, 0.75),
              .groups = "drop")
}

# for probability p, number of trials n, determine [q_low, q_up],
# so that a Bin(n, p) random variable is with probability nom_level in this interval
get_consistency_interval <- function(p, n, nom_level) {
  q_low <- qbinom((1 - nom_level) / 2, n, p, lower.tail = TRUE)
  q_up  <- qbinom((1 - nom_level) / 2, n, p, lower.tail = FALSE)
  data.frame(q_low, q_up) %>% 
    set_names(paste0(c("lower", "upper"), nom_level * 100))
}

# compute consistency bands for each quantile level
get_consistency_bands <- function(df) {
  df %>%
    group_by(model, quantile) %>%
    summarize(count = n(), .groups = "drop") %>%
    mutate(get_consistency_interval(quantile, count, 0.9) / count,
           get_consistency_interval(quantile, count, 0.5) / count) %>%
    select(-count)
}

plot_coverage <- function(df, 
                          date_column = target_end_date, 
                          B = 100, 
                          type = "confidence", 
                          difference = FALSE) {
  
  # labels
  labs <- c(consistency = "Consistency", confidence = "Confidence", confidence2 = "Confidence",
            Coverage = "Coverage", CoverageDiff = "Coverage - Level")
  # some customizations used in all plots
  my_theme <- list(
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = function(x) ifelse(x == 0, "0", x)),
    scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)),
    xlab("Quantile level"),
    theme_bw(base_size = 11),
    theme(panel.grid.major = element_line(size = 0.05), 
          panel.grid.minor = element_line(size = 0.05))
  )
  
  # compute coverage on full sample
  coverage_full <- coverage2(df)
  
  if (type %in% c("confidence", "confidence2")) {
    date_column <- enquo(date_column)
    bands <- get_confidence_bands(df, date_column, B)
    
    results <- coverage_full %>% 
      left_join(bands, by = c("model", "quantile"))
    
    if (!difference) {
      g <- ggplot(results) +
        facet_grid(type ~ model, labeller = labeller(type = labs)) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70") +
        {if (type == "confidence") list(
          geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkblue", alpha = 0.2),
          geom_ribbon(aes(x = quantile, ymin = u_5, ymax = u_95), fill = "darkred", alpha = 0.2))} +
        {if (type == "confidence2") list(
          geom_ribbon(aes(x = quantile, ymin = l_25, ymax = l_75), fill = "darkred", alpha = 0.3),
          geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkred", alpha = 0.2))} +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        ylab("Coverage") +
        coord_fixed()
    } else {
      results <- results %>% 
        mutate_at(vars("l", "u", "l_5", "l_25", "l_75", "l_95", "u_5", "u_25","u_75", "u_95"), 
                  list(~ . - quantile))
      
      g <- ggplot(results) +
        facet_grid(type ~ model, labeller = labeller(type = labs)) +
        geom_hline(yintercept = 0, size = 0.3, linetype = "solid", color = "darkgray") +
        {if (type == "confidence") list(
          geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkblue", alpha = 0.2),
          geom_ribbon(aes(x = quantile, ymin = u_5, ymax = u_95), fill = "darkred", alpha = 0.2))} +
        {if (type == "confidence2") list(
          geom_ribbon(aes(x = quantile, ymin = l_25, ymax = l_75), fill = "darkred", alpha = 0.3),
          geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkred", alpha = 0.2))} +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        ylab("Coverage - Level") +
        theme(aspect.ratio = 1)
    }
  } else if (type == "consistency") {
    bands <- get_consistency_bands(df)
    
    results <- coverage_full %>% left_join(bands, by = c("model", "quantile"))
    
    if (!difference) {
      g <- ggplot(results) +
        facet_grid(type ~ model, labeller = labeller(type = labs)) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        ylab("Coverage") +
        coord_fixed()
    } else {
      results <- results %>% 
        mutate_at(vars("l", "u", "lower50", "upper50", "lower90", "upper90"), 
                  list(~ . - quantile))
      
      g <- ggplot(results) +
        facet_grid(type ~ model, labeller = labeller(type = labs)) +
        geom_hline(yintercept = 0, size = 0.3, linetype = "solid", color = "darkgray") +
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        ylab("Coverage - Level") +
        theme(aspect.ratio = 1)
    }
  }
  # print(g)
  # invisible(results)
  return(g)
}

p1 <- plot_coverage(df1, B=100, type="consistency", difference=FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        axis.title.y = element_text(hjust = -0.35))
p2 <- plot_coverage(df2, B=100, type="confidence2", difference=FALSE) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        plot.background = element_blank()) +
  labs(y = element_blank())

p1 / p2

ggsave("figures/national_coverage.pdf", width=160, height=110, unit="mm", device = "pdf", dpi=300)


models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

df1 <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name == "Vermont",
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

df2 <- df1

df1$type <- "consistency"
df2$type <- "confidence"

p1 <- plot_coverage(df1, B=100, type="consistency", difference=FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text = element_text(size = 8),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        axis.title.y = element_text(hjust = -0.35))
p2 <- plot_coverage(df2, B=100, type="confidence", difference=FALSE) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        plot.background = element_blank()) +
  labs(y = element_blank())

p1 / p2

ggsave("figures/Vermont_coverage.pdf", width=160, height=110, unit="mm", device = "pdf", dpi=300)

### State level

models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

df1 <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name != "US",
         target == "1 wk ahead inc death",
         model %in% models) %>% 
  mutate(value = floor(value))

df2 <- df1

df1$type <- "Coverage"
df2$type <- "CoverageDiff"

p1 <- plot_coverage(df1, B=100, type="confidence", difference=FALSE) +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        strip.text.x = element_text(size = 8),
        #strip.text = element_text(size = 8),
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        #axis.title.y = element_text(hjust = -0.35)
        ) +
  labs(y = "Coverage")
p2 <- plot_coverage(df2, B=100, type="confidence", difference=TRUE) +
  theme(strip.background.x = element_blank(),
        strip.text.x = element_blank(),
        strip.text.y = element_blank(),
        plot.margin = unit(c(0, 0, 0, 0), "null"),
        plot.background = element_blank()) +
  labs(y = "Coverage - Level")
  #labs(y = element_blank())

p1 / p2

ggsave("figures/states_coverage_2.pdf", width=160, height=110, unit="mm", device = "pdf", dpi=300)
