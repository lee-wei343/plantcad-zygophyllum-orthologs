library(zoo)
library(tidyverse)

setwd("~/plantcad-zygophyllum-orthologs/")

msic <- read_csv("sandbox/msic-test2/out/test2_scores_table0.csv")
msic$pos <- seq_len(nrow(msic))
msic <- msic |>
  mutate(class = case_when(
    is.na(MSIC) ~ NA_character_,
    MSIC > 0.5 ~ "high",
    MSIC < -0.5 ~ "low",
    .default = "neutral"
  )) |>
  mutate(class = factor(class, levels = c("high", "neutral", "low")))
rolling <- rollapply(data = msic$MSIC**2,
                     width = list(-30:29),
                     FUN = mean,
                     partial = TRUE)
msic$mod_rolling <- (rolling) 

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
class_colors <- c("low" = palette[["blue"]],
                  "neutral" = palette[["gray"]],
                  "high" = palette[["vermillion"]])

ggplot(data = msic, aes(x = pos, y = MSIC)) +
  geom_col(aes(fill = class), width = 1L) +
  geom_hline(yintercept = 0) +
  scale_fill_manual(values = class_colors) +
  geom_line(aes(y = mod_rolling), linewidth = 0.3) +
  coord_cartesian(ylim = c(-1, 1)) +
  labs(y = "MSIC", x = "position") +
  theme_minimal()
