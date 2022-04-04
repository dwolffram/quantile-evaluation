library(tidyverse)
library(geomtextpath)
library(ggrepel)
source("reliability_functions.R")

# Load data
df <- read_csv("data/2022-01-03_df_processed.csv.gz", col_types = cols()) %>%
  filter(location == "US")

df <- df %>%
  filter(quantile %in% c(0.25, 0.5, 0.75),
         target == "1 wk ahead inc death",
         !model %in% c("USC-SI_kJalpha", "COVIDhub-4_week_ensemble", "COVIDhub_CDC-ensemble"))

results <- df %>%
  group_by(model, quantile) %>%
  do(reldiag(.$value, .$truth, alpha = unique(.$quantile), n_resamples = 99))

# summarize scores
scores <- results %>%
  group_by(quantile, model) %>%
  distinct(across(score:pval_ucond))

write.csv(scores, "data/2022-01-03_score_components_national.csv", row.names=FALSE)

scores <- read.csv("data/2022-01-03_score_components_states.csv")

# adjust model names (to save space)
scores$quantile <- as.factor(scores$quantile)
scores$model <- str_replace(scores$model, "COVIDhub-baseline", "Baseline")
scores$model <- as.character(lapply(strsplit(as.character(scores$model), "-"), '[[', 1))
scores$model <- str_replace(scores$model, "COVIDhub", "COVIDhub-ensemble")
scores$model <- as.factor(scores$model)
scores$model <- fct_relevel(scores$model, "Baseline", "COVIDhub-ensemble", "KITmetricslab")

# define isolines
iso <- scores %>%
  group_by(quantile) %>%
  summarize(intercept = seq(ceiling(max(dsc)) + first(unc)%%1 - ceiling(min(mcb)), # add decimal part of unc to ensure integer valued scores on isolines
                            -(ceiling(max(mcb))+ first(unc)%%1 - ceiling(min(mcb))), 
                            length.out = 20),
                            # by = -round((max(dsc) - min(dsc))/6)), 
            slope = 1, 
            unc = unique(unc),
            .groups = "drop") %>%
  mutate(score = round(unc - intercept), label = score)

# manually remove scores from isolines if there is overlap
# iso$label[c(1, 2, 7, 10, 25, 28, 32, 54, 55, 56, 66)] <- NA
iso$label[c(1, 2, 6, 8, 21, 23, 41, 42, 44)] <- NA

ggplot(data = scores) +
  facet_wrap('quantile', scales = "free", ncol = 3) +
  geom_abline(data = iso, aes(intercept = intercept, slope = slope), color = "lightgray", alpha = 0.5,
              size = 0.5) +
  geom_labelabline(data = iso, aes(intercept = intercept, slope = slope, label = label), color = "gray50",
                   hjust = 0.85, size = 7*0.36, text_only = TRUE, boxcolour = NA, straight = TRUE) +
  geom_point(aes(x = mcb, y = dsc, color = model), size = 0.4) +
  geom_text_repel(aes(x = mcb, y = dsc, label = model), max.overlaps = NA, size = 8*0.36, nudge_x = 0.1,
                  direction = "both", segment.color = "transparent", box.padding = 0.1, force = 1, point.padding = 0) +
  xlab("MCB") +
  ylab("DSC") +
  #labs(title = "1 wk ahead inc case") +
  theme_bw(base_size = 11) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        aspect.ratio = 1,
        legend.position = "none",
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank()
  ) + 
  scale_color_brewer(palette="Set1")

ggsave("figures/6_score_decomposition_states.pdf", width=160, height=70, unit="mm", device = "pdf", dpi=300)
