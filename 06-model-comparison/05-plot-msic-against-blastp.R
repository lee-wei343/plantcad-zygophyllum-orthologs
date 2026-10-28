library(tidyverse)
library(ggplot2)


per_tx <- msic_blastp |>
  filter(in_cds) |>
  filter(!model %in% c("Botanic1L", "Botanic1S", "PCAD2L")) |>
  summarise(
    mean_msic = mean(msic, na.rm = TRUE),
    pident = first(pident),
    .by = c(species, group, model, tx)
  )

ggplot(data = per_tx, aes(x = pident, y = mean_msic, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(group ~ model) +
  theme_minimal()
