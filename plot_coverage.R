library(tidyverse)
Sys.setlocale("LC_ALL", "C")

coverage <- function(df) {
  df$l <- df$truth < floor(df$value)
  df$u <- df$truth <= floor(df$value)
  
  df %>%
    group_by(model, quantile) %>%
    summarize(l = mean(l), u = mean(u), .groups = "drop")
}

# get critical value of binomial test (used for consistency bands)
get_cr <- function(sample_size, alpha, nominal_level, alternative = "less") {
  if (alternative == "less") {
    C <- 0
    while (pbinom(C + 1, sample_size, alpha) < nominal_level) C <- C + 1
  } else if (alternative == "greater") {
    C <- sample_size
    while (1 - pbinom(C - 1, sample_size, alpha) < nominal_level) C <- C - 1
  }
  return(C)
}

sample_n_groups = function(grouped_df, n, replace = FALSE) {
  groups_resampled <- grouped_df %>% 
    group_keys() %>% 
    slice_sample(n = n, replace = replace)
  
  grouped_df %>% 
    right_join(groups_resampled, by = group_vars(grouped_df))
}

plot_coverage <- function(df, 
                          date_column = target_end_date, 
                          B = 1000, 
                          type = "confidence", 
                          difference = FALSE) {
  
  # some customizations used in all plots
  my_theme <- list(
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = function(x) ifelse(x == 0, "0", x)),
    scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)),
    xlab("Quantile level"),
    theme_bw(base_size = 11),
    theme(panel.grid.major = element_line(size = 0.05), 
          panel.grid.minor = element_line(size = 0.05))
  )
  
  if (type %in% c("confidence", "confidence2")) {
    date_column <- enquo(date_column)

    # compute coverage on all bootstrap samples
    coverage_df = data.frame() 
    
    for (i in 1:B) {
      df_resampled <- df %>%
        group_by(!!date_column) %>%
        sample_n_groups(n = n_distinct(group_keys(.)), replace = TRUE)
      
      coverage_sample <- coverage(df_resampled)
      
      coverage_df <- bind_rows(coverage_df, coverage_sample)
    }
    
    # compute CIs from bootstrapped coverage
    results <- coverage_df %>%
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
    
    # compute coverage on full sample
    coverage_full <- coverage(df)
    
    results <- results %>%
      left_join(coverage_full, by = c("model", "quantile"))
    
    if (!difference) {
      g <- ggplot(results) +
        facet_wrap("model", ncol = 3) +
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
    }
    
    else {
      results <- results %>% 
        mutate_at(vars("l", "u", 
                       "l_5", "l_25", "l_75", "l_95", 
                       "u_5", "u_25","u_75", "u_95"), 
                  list(~ . - quantile))
      
      g <- ggplot(results) +
        facet_wrap("model") +
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
  }
  
  else if (type == "consistency") {
    results <- coverage(df)
    
    consistency_bands <- df %>% 
      group_by(model, quantile) %>% 
      summarize(count = n(), 
                .groups = "drop")
    
    consistency_bands <- consistency_bands %>% 
      rowwise() %>% 
      mutate(lower90 = get_cr(count, quantile, 0.05, "less")/count,
             upper90 = get_cr(count, quantile, 0.05, "greater")/count,
             lower50 = get_cr(count, quantile, 0.25, "less")/count,
             upper50 = get_cr(count, quantile, 0.25, "greater")/count)
    
    results <- results %>% 
      left_join(consistency_bands, by = c("model", "quantile"))
    
    if (!difference) {
      g <- ggplot(results) +
        facet_wrap("model", ncol = 3) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        ylab("Coverage") +
        coord_fixed()
    }
    
    else {
      results <- results %>% 
        mutate_at(vars("l", "u", "lower50", "upper50", "lower90", "upper90"), 
                  list(~ . - quantile))
      
      g <- ggplot(results) +
        facet_wrap("model") +
        geom_hline(yintercept = 0, size = 0.3, linetype = "solid", color = "darkgray") +
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        ylab("Coverage - Level") +
        theme(aspect.ratio = 1)
    }
  }
  print(g)
  invisible(results)
}



# Load data
models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

# National level

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location == "US")

df <- df %>%
  filter(target == "1 wk ahead inc death",
         model %in% models)

plot_coverage(df, B=100, type="confidence2", difference=FALSE)
ggsave("figures/national_coverage_confidence.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

plot_coverage(df, B=100, type="consistency", difference=FALSE)
ggsave("figures/national_coverage_consistency.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

# State level

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name != "US")

df <- df %>%
  filter(target == "1 wk ahead inc death",
         model %in% models)

plot_coverage(df, B=100, type="confidence", difference=FALSE)
ggsave("figures/states_coverage.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

plot_coverage(df, B=100, type="confidence", difference=TRUE)
ggsave("figures/states_coverage_diff.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)


# Vermont

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name == "Vermont")

df <- df %>%
  filter(target == "1 wk ahead inc death",
         model %in% models)

plot_coverage(df, B=100, type="confidence", difference=FALSE)
ggsave("figures/Vermont_coverage_confidence.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)

plot_coverage(df, B=100, type="consistency", difference=FALSE)
ggsave("figures/Vermont_coverage_consistency.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)






### Multiple states at once

models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name %in% c("Nebraska", "Vermont", "Florida"))

df <- df %>%
  filter(target == "1 wk ahead inc death",
         model %in% models)

coverage2 <- function(df) {
  df$l <- df$truth < floor(df$value)
  df$u <- df$truth <= floor(df$value)
  
  df %>%
    group_by(model, location_name, quantile) %>%
    summarize(l = mean(l), u=mean(u), .groups = "drop")
}

coverage_df <- coverage2(df)

consistency_bands <- df %>% 
  group_by(model, location_name, quantile) %>% 
  summarize(count = n())

consistency_bands <- consistency_bands %>% 
  rowwise() %>% 
  mutate(lower90 = get_cr(count, quantile, 0.05, "less")/count,
         upper90 = get_cr(count, quantile, 0.05, "greater")/count,
         lower50 = get_cr(count, quantile, 0.25, "less")/count,
         upper50 = get_cr(count, quantile, 0.25, "greater")/count)

coverage_df <- coverage_df %>% 
  left_join(consistency_bands)

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
