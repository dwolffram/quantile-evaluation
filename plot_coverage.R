library(tidyverse)
Sys.setlocale("LC_ALL", "C")

coverage <- function(df){
  df$l <- df$truth < floor(df$value)
  df$u <- df$truth <= floor(df$value)
  
  df <- df %>%
    group_by(model, quantile) %>%
    summarize(l = mean(l), u=mean(u), .groups = "drop")
  
  return(df)
}

# get critical value of binomial test (used for consistency bands)
get_cr <- function(sample_size, alpha, nominal_level, alternative='less'){
  if(alternative == 'less'){
    C <- 0
    while(pbinom(C + 1, sample_size, alpha) < nominal_level) C <- C+1
  } 
  else if(alternative == 'greater'){
    C <- sample_size
    while (1 - pbinom(C-1, sample_size, alpha) < nominal_level) C <- C -1
  }
  return(C)
}

plot_coverage <- function(df, date_column = target_end_date, B = 1000, type = "confidence", difference = FALSE){
  
  # some customizations used in all plots
  my_theme <- list(
    scale_x_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1), labels = function(x) ifelse(x == 0, "0", x)),
    scale_y_continuous(labels = function(y) ifelse(y == 0, "0", y)),
    xlab("Quantile Level"),
    ylab(NULL),
    theme_bw(base_size = 11),
    theme(panel.grid.major = element_line(size = 0.05), 
          panel.grid.minor = element_line(size = 0.05))
  )
  
  if (type == "confidence"){
    date_column <- enquo(date_column)
    dates <- df %>% distinct(!!date_column) %>% pull()
    
    # compute coverage on all bootstrap samples
    coverage_df = data.frame()  
    for(i in 1:B){
      dates_resampled <- sample(dates, replace = TRUE)
      
      # simple filter doesn't work because the same date can occur multiple times
      df_resampled <- lapply(dates_resampled, function(x){df %>% filter(target_end_date == x)}) %>%
        bind_rows()
      
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
    
    if (!difference){
      g <- ggplot(results) +
        facet_wrap("model", ncol = 3) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
        geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkblue", alpha = 0.2) +
        geom_ribbon(aes(x = quantile, ymin = u_5, ymax = u_95), fill = "darkred", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
        coord_fixed()
    }
    
    else{
      results <- results %>% 
        mutate_at(vars("l", "u", 
                       "l_5", "l_25", "l_75", "l_95", 
                       "u_5", "u_25","u_75", "u_95"), 
                  list(~ . - quantile))
      
      g <- ggplot(results) +
        facet_wrap("model") +
        geom_hline(yintercept = 0, size = 0.3, linetype = "solid", color = "darkgray") +
        geom_ribbon(aes(x = quantile, ymin = l_5, ymax = l_95), fill = "darkblue", alpha = 0.2) +
        geom_ribbon(aes(x = quantile, ymin = u_5, ymax = u_95), fill = "darkred", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme
    }
  }
  
  else if (type == "consistency"){
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
    
    if (!difference){
      g <- ggplot(results) +
        facet_wrap("model", ncol = 3) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 0.2, linetype = "solid", colour = "grey70")+
        geom_ribbon(aes(x = quantile, ymin = lower50, ymax = upper50), fill = "skyblue3", alpha = 0.3) +
        geom_ribbon(aes(x = quantile, ymin = lower90, ymax = upper90), fill = "skyblue3", alpha = 0.2) +
        geom_errorbar(aes(x=quantile, ymin = l, ymax = u), width = 0.0125, size = 0.3, colour = "black") +
        my_theme +
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
        my_theme
      
    }
  }
  print(g)
  invisible(results)
}



# Load data
models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

# df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
#   filter(location == "US")

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name == "Vermont")

df <- df %>%
  filter(target == "1 wk ahead inc death",
         model %in% models)

# Plot coverage
results <- plot_coverage(df, B = 100)
results_diff <- plot_coverage(df, B = 100, difference = TRUE)

plot_coverage(df, B=100, type="confidence", difference=FALSE)
plot_coverage(df, B=100, type="confidence", difference=TRUE)
plot_coverage(df, B=100, type="consistency", difference=FALSE)
plot_coverage(df, B=100, type="consistency", difference=TRUE)


# ggsave("figures/coverage_national_consistency.pdf", width=180, height=70, unit="mm", device = "pdf", dpi=300)




### Multiple states at once

models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location_name %in% c("Nebraska", "Vermont", "Florida"))

df <- df %>%
  filter(target == "1 wk ahead inc death",
         model %in% models)

coverage2 <- function(df){
  df$l <- df$truth < floor(df$value)
  df$u <- df$truth <= floor(df$value)
  
  df <- df %>%
    group_by(model, location_name, quantile) %>%
    summarize(l = mean(l), u=mean(u), .groups = "drop")
  
  return(df)
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
  xlab("Quantile") +
  ylab(NULL) +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_line(size = 0.05)) +
  coord_fixed()

# ggsave("figures/coverage_3states_consistency.pdf", width=200, height=200, unit="mm", device = "pdf", dpi=300)
