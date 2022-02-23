library(tidyverse)

Sys.setlocale("LC_ALL", "C")

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols())
df_wide <- pivot_wider(df, names_from=quantile, names_prefix="value.", values_from=value)

df_ens <- df %>%
  filter(model == "COVIDhub-ensemble",
         str_detect(target, "death"),
         location == "US")

cols <- colorRampPalette(c("deepskyblue4", "lightgrey"))(2 + 1)[-1]

ggplot(df_ens, aes(x=target_end_date, y=truth)) +
  facet_wrap("target", scales="free_y", ncol=2) +
  geom_smooth(aes(y = value.0.5, ymin = value.0.025, ymax = value.0.975), 
              linetype=3, size=0.7, colour="white", fill=cols[2], alpha=0.8, stat = "identity") +
  geom_smooth(aes(y = value.0.5, ymin = value.0.25, ymax = value.0.75),
              linetype=3, size=0.7, colour="white", fill=cols[1], alpha=0.8, stat = "identity") +
  geom_line() +
  geom_point(pch = 4, size=3) +
  geom_point(aes(y = value.0.5), pch = 21, col = "black", bg = "white", size=3) +
  theme_bw() +
  labs(title="National level forecasts - KITmetricslab-select_ensemble",
       x = NULL,
       y = "Deaths") +
  theme_bw(base_size=16) +
  theme(strip.text.y = element_text(size = 8))



df2 <- df %>%
  filter(model == "COVIDhub-ensemble",
         target == "1 wk ahead inc death",
         location == "US",
         !quantile %in% c(0.01, 0.5, 0.99) )

df2_wide <- df_wide %>%
  filter(model == "COVIDhub-ensemble",
         target == "1 wk ahead inc death",
         location == "US")

ggplot(df2_wide, aes(x=target_end_date)) +
  geom_point(data=df2, aes(x=target_end_date, y=value), 
             shape='-', size=3, col='deepskyblue3', alpha=0.9) +
  geom_crossbar(aes(y=value.0.5, ymin=value.0.01, ymax=value.0.99), fatten=1.5,
                width=3, size=0.5, col='deepskyblue4', fill='gray', alpha=0.3) +
  geom_line(aes(y=truth), size=1, col='darkred') +
  scale_x_date(date_breaks = "months" , date_labels = "%b-%Y") +
  xlab(NULL) +
  ylab('Incident deaths') +
  labs(title = 'One-week-ahead ensemble forecasts') +
  theme_bw(base_size = 11)

ggsave("figures/ensemble_forecast.pdf", width=180, height=100, unit="mm", device = "pdf", dpi=300)


### Multiple models

unique(df$model)

models <- c('COVIDhub-baseline', 'COVIDhub-ensemble', 'KITmetricslab-select_ensemble')

df2 <- df %>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US",
         !quantile %in% c(0.01, 0.5, 0.99) )

df2_wide <- df_wide %>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US")

ggplot(df2_wide, aes(x=target_end_date)) +
  facet_wrap("model", ncol = 1) +
  geom_point(data=df2, aes(x=target_end_date, y=value), 
             shape='-', size=1.8, col='darkslategray4', alpha=0.9) +
  geom_crossbar(aes(y=value.0.5, ymin=value.0.01, ymax=value.0.99), fatten=4,
                width=2, size=0.1, col='deepskyblue4', fill='gray', alpha=0.3) +
  geom_line(aes(y=truth), size=0.4, col='darkred') +
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y") +
  xlab(NULL) +
  ylab('Incident deaths') +
  labs(title = 'One-week-ahead forecasts') +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_line(size = 0.05))

ggsave("figures/ensemble_forecast3.pdf", width=180, height=180, unit="mm", device = "pdf", dpi=500)


#### 

highlight <- c(0.05, 0.25, 0.50, 0.75, 0.95)

df2 <- df %>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US",
         !quantile %in% c(0.01, 0.5, 0.99) )

df3 <- df2 %>%
  filter(quantile %in% highlight)

df2 <- df2 %>%
  filter(!quantile %in% highlight)

df2_wide <- df_wide %>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US")

ggplot(df2_wide, aes(x=target_end_date)) +
  facet_wrap("model", ncol = 1) +
  geom_point(data=df2, aes(x=target_end_date, y=value), 
             shape='-', size=0.5, col='black', alpha=0.9) +
  # geom_point(data=df3, aes(x=target_end_date, y=value), 
  #            shape='-', size=2, col='deepskyblue4') +
  geom_crossbar(aes(y=value.0.5, ymin=value.0.01, ymax=value.0.99), fatten=4,
                width=2, size=0.1, col='black', fill='deepskyblue4', alpha=0.1) +
  geom_crossbar(aes(y=value.0.5, ymin=value.0.05, ymax=value.0.95), fatten=4,
                width=2, size=0.1, col='black', fill='deepskyblue4', alpha=0.4) +
  geom_crossbar(aes(y=value.0.5, ymin=value.0.25, ymax=value.0.75), fatten=4,
                width=2, size=0.1, col='black', fill='deepskyblue4', alpha=0.8) +
  geom_line(aes(y=truth), size=0.4, col='darkred') +
  scale_x_date(date_breaks = "months" , date_labels = "%b-%y") +
  xlab(NULL) +
  ylab('Incident deaths') +
  labs(title = 'One-week-ahead forecasts') +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_line(size = 0.05))

ggsave("figures/ensemble_forecast3.pdf", width=180, height=180, unit="mm", device = "pdf", dpi=500)


###

models <- c('COVIDhub-baseline', 'COVIDhub-ensemble', 'KITmetricslab-select_ensemble')

df2 <- df %>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US"#,
         #!quantile %in% c(0.01, 0.5, 0.99) 
         )

df2_wide <- df_wide %>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US")

highlight <- c(0.05, 0.25, 0.50, 0.75, 0.95)

df2 <- df2 %>% 
  mutate(color = quantile %in% highlight)

ggplot(df2_wide, aes(x=target_end_date)) +
  facet_wrap("model", ncol = 1) +
  geom_crossbar(aes(y=value.0.5, ymin=value.0.01, ymax=value.0.99), fatten=1,
                width=2, size = 0.3, colour = "gray") +
  geom_segment(data = df2, aes(x = target_end_date - 1, xend = target_end_date + 1, y = value, yend = value,
                               colour = color), size = 0.3) + 
  geom_line(aes(y=truth), size = 0.3, col='darkred') +
  scale_x_date(date_breaks = "months" , date_labels = "%b", expand = c(0.02, 0)
              # limits = c(as.Date("2021-04-20"), as.Date("2022-01-5"))
               ) +
  xlab(NULL) +
  ylab('Incident deaths') +
  # labs(title = 'One-week-ahead forecasts') +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(hjust = -1)) +
  scale_color_manual(values = c("gray", "black"))

ggsave("figures/forecasts.pdf", width=200, height=220, unit="mm", device = "pdf", dpi=500)

highlight <- c(0.025, 0.25, 0.50, 0.75, 0.975)

df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols())%>%
  filter(model %in% models,
         target == "1 wk ahead inc death",
         location == "US")

df1 <- df %>% 
  filter(quantile %in% highlight)

df2 <- df %>%
  filter(!quantile %in% highlight)

df_box <- df %>%
  filter(quantile %in% c(0.01, 0.5, 0.99)) %>% 
  pivot_wider(names_from = quantile, names_prefix = "value.", values_from = value)

ggplot(df1, aes(x=target_end_date)) +
  facet_wrap("model", ncol = 1) +
  geom_crossbar(data = df_box, aes(y = value.0.5, ymin = value.0.01, ymax = value.0.99), fatten = 1,
                width = 2, size = 0.3, colour = "gray") +
  geom_segment(data = df2, aes(x = target_end_date - 1, xend = target_end_date + 1, y = value, yend = value), 
               color = "gray", size = 0.3) + 
  geom_segment(data = df1, aes(x = target_end_date - 1, xend = target_end_date + 1, y = value, yend = value), 
               color = "black", size = 0.4) + 
  geom_line(aes(y=truth), size = 0.3, col='darkred') +
  scale_x_date(date_breaks = "months" , date_labels = "%b", expand = c(0.02, 0)) +
  xlab(NULL) +
  ylab('Incident deaths') +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(hjust = -1.25))

ggsave("figures/forecasts.pdf", width=180, height=220, unit="mm", device = "pdf", dpi=500)

ggplot(df1, aes(x=target_end_date)) +
  facet_wrap("model", ncol = 1) +
  geom_crossbar(data = df_box, aes(y = value.0.5, ymin = value.0.01, ymax = value.0.99), fatten = 1,
                width = 2, size = 0.3, colour = "gray", fill = "gray", alpha = 0.3) +
  geom_segment(data = df2, aes(x = target_end_date - 1, xend = target_end_date + 1, y = value, yend = value), 
               color = "gray", size = 0.3) + 
  geom_segment(data = df1, aes(x = target_end_date - 1.25, xend = target_end_date + 1.25, y = value, yend = value), 
               color = "deepskyblue4", size = 0.4, lineend = "round") + 
  geom_line(aes(y=truth), size = 0.3, col='darkred') +
  scale_x_date(date_breaks = "months" , date_labels = "%b", expand = c(0.02, 0)) +
  xlab(NULL) +
  ylab('Incident deaths') +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_line(size = 0.05), 
        panel.grid.minor = element_blank(),
        legend.position = "none",
        axis.text.x = element_text(hjust = -1.25))

ggsave("figures/forecasts2.pdf", width=180, height=220, unit="mm", device = "pdf", dpi=500)

