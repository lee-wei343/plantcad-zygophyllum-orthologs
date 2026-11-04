library(tidyverse)
library(ggplot2)

setwd("~/plantcad-zygophyllum-orthologs")

plot_dir <- file.path("06-model-comparison", "out")
dir.create(plot_dir, showWarnings = FALSE)

msic_blastp <- read_csv("06-model-comparison/msic_blastp.csv",
                        show_col_types = FALSE)

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
species_colors <- c(
  "a_thaliana" = palette[["blue"]],
  "z_stapffii" = palette[["vermillion"]]
)

# -------------------------------------------------------------------------


per_tx <- msic_blastp |>
  filter(in_cds) |>
  filter(!model %in% c("Botanic1L", "Botanic1S", "PCAD2L", "PCAD2S")) |>
  summarise(
    mean_msic = mean(msic, na.rm = TRUE),
    pident = first(pident),
    .by = c(species, group, model, pair, tx)
  )

p1 <- ggplot(data = per_tx, aes(x = pident, y = mean_msic, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  scale_color_manual(values = species_colors) +
  labs(y = "Average MSIC per transcript",
       x = "Percentage identity of best BLASTp segment per transcript")

ggsave(file.path(plot_dir, "msic_x_blastp_pcad1_exons.png"),
       plot = p1,
       dpi = 300,
       height = 5,
       width = 8)


per_tx_all <- msic_blastp |>
  filter(in_cds) |>
  summarise(
    mean_msic = mean(msic, na.rm = TRUE),
    pident = first(pident),
    .by = c(species, group, model, pair, tx)
  )

all <- ggplot(data = per_tx_all, aes(x = pident, y = mean_msic, colour = species)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  facet_wrap(model~group) +
  scale_color_manual(values = species_colors) +
  theme_minimal() +
  labs(y = "Average MSIC per transcript",
       x = "Percentage identity of best BLASTp segment per transcript")

ggsave(file.path(plot_dir, "msic_x_blastp_all_exons.png"),
       plot = all,
       dpi = 300,
       height = 5,
       width = 8)


# atMSIC_x_zsMSIX scatter -------------------------------------------------

pair <- per_tx |>
  select(pair, species, mean_msic) |>
  pivot_wider(names_from = species, values_from = mean_msic)

ggplot(data = pair, aes(x = a_thaliana, y = z_stapffii)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE) +
  theme_minimal() +
  labs(x = "A. thaliana", y = "Z. stapffii", title = "Mean MSIC per Transcript")

c <- cor.test(pair$z_stapffii, pair$a_thaliana, method = "spearman")

