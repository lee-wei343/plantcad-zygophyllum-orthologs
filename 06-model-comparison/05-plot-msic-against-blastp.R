library(tidyverse)
library(ggplot2)

setwd("~/plantcad-zygophyllum-orthologs")

msic_blastp <- read_csv("06-model-comparison/msic_blastp.csv",
                        show_col_types = FALSE)

per_tx <- msic_blastp |>
  filter(in_cds) |>
  filter(!model %in% c("Botanic1L", "Botanic1S", "PCAD2L", "PCAD2S")) |>
  summarise(
    mean_msic = mean(msic, na.rm = TRUE),
    pident = first(pident),
    .by = c(species, group, model, pair, tx)
  )

ggplot(data = per_tx, aes(x = pident, y = mean_msic, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal()


pair <- per_tx |>
  select(pair, species, mean_msic) |>
  pivot_wider(names_from = species, values_from = mean_msic)

ggplot(data = pair, aes(x = a_thaliana, y = z_stapffii)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(x = "A. thaliana", y = "Z. stapffii", title = "Mean MSIC per Transcript")

c <- cor.test(pair$z_stapffii, pair$a_thaliana, method = "spearman")

c
