library(tidyverse)
library(gridExtra)
Sys.setlocale("LC_ALL", "C")

elementary_quantile_score <- function(y_true, y_pred, theta, alpha){
  ((y_true < y_pred) - alpha) * ((theta < y_pred) - (theta < y_true))
}

get_thetas <- function(df, n=1001){
  tmp <- c(df$value, df$truth)
  thetas <- seq(from = min(tmp) - 1, to = max(tmp) + 1, length.out = n)
  return(thetas)
}

quantile_score <- function(y_true, y_pred, alpha){
  if (y_true < y_pred){
    (1 - alpha)*(y_pred - y_true)
  }
  else{
    alpha*(y_true - y_pred)
  }
}

murphy_diagram <- function(df){
  scores <- df %>% 
    rowwise() %>% 
    mutate(qs = quantile_score(truth, value, quantile)) %>% 
    group_by(model, quantile) %>% 
    summarize(mean_qs = mean(qs), 
              .groups = "drop")
  
  
  scores <- scores %>% 
    group_by(quantile, model) %>% 
    summarize(label = paste0(model, " (", round(mean_qs, digits = 1), ")", collapse = "\n"), 
              .groups = "drop")
  
  
  df <- df %>% 
    left_join(df %>%
                group_by(quantile) %>%
                do(crossing(theta = get_thetas(., 150))), by = "quantile")
  
  df <- df %>%
    mutate(score = elementary_quantile_score(truth, value, theta, quantile))
  
  df <- df %>%
    group_by(model, quantile, theta) %>%
    summarize(mean_score = mean(score), .groups="keep")
  
  df <- df %>% 
    left_join(scores, by = c("model", "quantile"))
  
  
  ymax = max(df$mean_score)
  xmax = max(df$theta)
  
  xs <- split(df,f = df$quantile)
  p1 <- ggplot(xs$`0.25`) +
    geom_line(aes(x=theta, y=mean_score, color=label), size=0.5) +
    facet_wrap("quantile", scales = "free") +
    xlab(expression(paste("Threshold ", theta))) +
    ylab("Elementary score") +
    theme_bw(base_size = 11) +
    theme(aspect.ratio = 1, 
          legend.justification=c(1,1), legend.position=c(0.95,1), 
          legend.title=element_text(size=6, face = "bold"),
          legend.text=element_text(size=6),
          legend.title.align = 0, 
          legend.text.align = 0,
          legend.key.size = unit(0.2, "lines"),
          legend.background = element_blank(),
          panel.grid.major = element_line(size = 0.05), 
          panel.grid.minor = element_line(size = 0.05)) + 
    scale_color_brewer(palette="Set1") +
    labs(color = "Model (quantile score)") +
    expand_limits(x = xmax, y = ymax)
  
  p2 <- p1 %+% xs$`0.5`
  p3 <- p1 %+% xs$`0.75`
  
  g <- arrangeGrob(p1 + xlab(""), p2, p3 + xlab(""), ncol = 3)
  plot(g)
  invisible(g)
}


# Load and filter data
df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location != "US")

models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")
target <- "1 wk ahead inc death"
quantiles <- c(0.25, 0.5, 0.75)

df <- df %>%
  filter(model %in% models,
         target == !!target,
         quantile %in% quantiles)

# Shorten model names
df$model <- str_replace(df$model, "COVIDhub-baseline", "Baseline")
df$model <- as.character(lapply(strsplit(as.character(df$model), "-"), '[[', 1))
df$model <- str_replace(df$model, "COVIDhub", "COVIDhub-ensemble")
df$model <- as.factor(df$model)
df$model <- fct_relevel(df$model, "Baseline", "COVIDhub-ensemble", "KITmetricslab")

# Plot murphy diagram
g <- murphy_diagram(df)

# ggsave("figures/murphy_states_qs.pdf", plot=g, width=180, height=70, unit="mm", device = "pdf", dpi=300)




### OLD CODE

# murphy_diagram <- function(df, models, target, 
#                            quantiles = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975, 0.01, 0.05, 0.15, 0.2, 
#                                          0.3, 0.35, 0.4, 0.45, 0.55, 0.6, 0.65, 0.7, 0.8, 0.85, 0.95, 0.99)){
#   df <- df %>%
#     filter(model %in% models,
#            target == !!target,
#            quantile %in% quantiles)
#   
#   df <- df %>% left_join(df %>%
#                            group_by(quantile) %>%
#                            do(crossing(theta = get_thetas(., 150))), by = "quantile")
#   
#   df <- df %>%
#     mutate(score = elementary_quantile_score(truth, value, theta, quantile))
#   
#   df <- df %>%
#     group_by(model, quantile, theta) %>%
#     summarize(mean_score = mean(score), .groups="keep")
#   
#   ggplot(df, aes(x=theta, y=mean_score, color=model)) +
#     geom_line(size=0.5) +
#     facet_wrap("quantile", scales = "free") +
#     # labs(title = paste("Murphy diagram:", target), color = "Model") +
#     xlab(expression(paste("Parameter ", theta))) +
#     ylab(NULL) +
#     theme_bw(base_size = 11) +
#     theme(aspect.ratio = 1, legend.position = "bottom", legend.title = element_blank(),
#           legend.background = element_blank(),
#           legend.box.background = element_rect(color = "black"),
#           panel.grid.major = element_line(size = 0.05), 
#           panel.grid.minor = element_line(size = 0.05)) + 
#     scale_color_brewer(palette="Set1")
# }

# df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
#   filter(location != "US")
# 
# models <- c("KITmetricslab-select_ensemble", "COVIDhub-ensemble", "COVIDhub-baseline")
# target <- "1 wk ahead inc death"
# quantiles <- c(0.25, 0.5, 0.75)
# 
# df <- df %>%
#   filter(model %in% models,
#          target == !!target,
#          quantile %in% quantiles)
# 
# 
# scores <- df %>% 
#   rowwise() %>% 
#   mutate(qs = quantile_score(truth, value, quantile)) %>% 
#   group_by(model, quantile) %>% 
#   summarize(mean_qs = mean(qs))
# 
# 
# scores <- scores %>% 
#   group_by(quantile) %>% 
#   summarize(label = paste0("Quantile score\n", paste0(model, ": ", round(mean_qs, digits = 1), collapse = "\n")))
# 
# 
# df <- df %>% left_join(df %>%
#                          group_by(quantile) %>%
#                          do(crossing(theta = get_thetas(., 150))), by = "quantile")
# 
# df <- df %>%
#   mutate(score = elementary_quantile_score(truth, value, theta, quantile))
# 
# df <- df %>%
#   group_by(model, quantile, theta) %>%
#   summarize(mean_score = mean(score), .groups="keep")
# 
# ggplot(df) +
#   geom_line(aes(x=theta, y=mean_score, color=model), size=0.5) +
#   facet_wrap("quantile", scales = "free") +
#   geom_label(data = scores, mapping = aes(x = Inf, y = Inf, label = label),
#              size = 6*0.36, hjust = 1, vjust = 1, label.size = NA, alpha=0, label.padding = unit(1, "lines")) +
#   # labs(title = paste("Murphy diagram:", target), color = "Model") +
#   xlab(expression(paste("Parameter ", theta))) +
#   ylab(NULL) +
#   theme_bw(base_size = 11) +
#   theme(aspect.ratio = 1, 
#         legend.position = "bottom", legend.title = element_blank(),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(color = "black"),
#         panel.grid.major = element_line(size = 0.05), 
#         panel.grid.minor = element_line(size = 0.05)) + 
#   scale_color_brewer(palette="Set1")

# ggsave("figures/murphy_states_qs.pdf", width=180, height=100, unit="mm", device = "pdf", dpi=300)
