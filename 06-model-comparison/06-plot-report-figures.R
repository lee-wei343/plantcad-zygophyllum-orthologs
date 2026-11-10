library(tidyverse)

# Report figures for the Results section

setwd("~/plantcad-zygophyllum-orthologs")

fig_dir <- file.path("report", "Figures")
dir.create(fig_dir, showWarnings = FALSE)

font <- "Libertinus Sans"
theme_set(theme_minimal(base_size = 9, base_family = font))

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
species_colors <- c(
  "A. thaliana" = palette[["blue"]],
  "Z. stapffii" = palette[["vermillion"]]
)
model_levels <- c("PCAD1", "PCAD2S", "PCAD2L", "Botanic1S", "Botanic1L")
model_labels <- c(
  PCAD1 = "PlantCAD1", PCAD2S = "PlantCAD2-S", PCAD2L = "PlantCAD2-L",
  Botanic1S = "Botanic1-S", Botanic1L = "Botanic1-L"
)
thr <- 0.5


# load --------------------------------------------------------------------

cds <- data.table::fread(
  "06-model-comparison/msic_long.csv",
  select = c("tx", "species", "model", "msic", "in_cds", "group")
) |>
  as_tibble() |>
  filter(in_cds, !is.na(msic)) |>
  mutate(
    species = recode(species, a_thaliana = "A. thaliana",
                     z_stapffii = "Z. stapffii"),
    gene_set = recode(group, top15 = "BUSCO", tps = "TPS"),
    model = factor(model, levels = model_levels, labels = model_labels)
  )

best_blastp <- function(dir_path) {
  list.files(dir_path, pattern = "_with_.*\\.tsv$") |>
    set_names() |>
    map(\(f) read_tsv(file.path(dir_path, f), show_col_types = FALSE)) |>
    list_rbind(names_to = "file") |>
    mutate(
      at_tx = str_extract(file, "^[^_]+"),
      zs_tx = str_extract(file, "_with_(.+)\\.tsv", group = 1)
    ) |>
    slice_max(bitscore, n = 1, with_ties = FALSE, by = c(at_tx, zs_tx)) |>
    select(at_tx, zs_tx, pident)
}

pairs <- bind_rows(
  best_blastp(file.path("04-top15", "out", "blastp")),
  best_blastp(file.path("05-additional-genes", "blastp"))
)

per_tx <- cds |>
  summarise(mean_msic = mean(msic), .by = c(species, gene_set, model, tx))

per_pair <- pairs |>
  left_join(
    per_tx |> filter(species == "A. thaliana") |>
      select(model, gene_set, at_tx = tx, at = mean_msic),
    by = "at_tx"
  ) |>
  left_join(
    per_tx |> filter(species == "Z. stapffii") |>
      select(model, zs_tx = tx, zs = mean_msic),
    by = c("zs_tx", "model")
  )


# fig: CDS MSIC distributions ---------------------------------------------

labels <- cds |>
  mutate(region = cut(msic, c(-Inf, -thr, thr, Inf),
                      labels = c("low", "mid", "high"))) |>
  count(model, species, region) |>
  mutate(
    label = scales::percent(n / sum(n), accuracy = 0.1),
    x = c(low = -thr, mid = 0, high = thr)[as.character(region)],
    hjust = c(low = 1.1, mid = 0.5, high = -0.1)[as.character(region)],
    .by = c(model, species)
  )

p_hist <- ggplot(cds, aes(x = msic, fill = species)) +
  geom_histogram(
    aes(y = after_stat(count / ave(count, PANEL, FUN = sum))),
    binwidth = 0.05, boundary = 0, show.legend = FALSE
  ) +
  geom_vline(xintercept = c(-thr, thr), linewidth = 0.3) +
  geom_text(data = labels,
            aes(x = x, y = Inf, label = label, hjust = hjust),
            vjust = 1.3, size = 2.4, family = font, inherit.aes = FALSE) +
  facet_grid(model ~ species) +
  scale_fill_manual(values = species_colors) +
  scale_y_continuous(labels = scales::label_percent(),
                     expand = expansion(mult = c(0, 0.35))) +
  labs(x = "MSIC", y = "Relative frequency of CDS bases") +
  theme(strip.text.x = element_text(face = "italic"))

ggsave(file.path(fig_dir, "msic_hist_cds.pdf"), p_hist,
       width = 16, height = 18, units = "cm", device = cairo_pdf)


paired_long <- per_pair |>
  pivot_longer(c(at, zs), names_to = "species", values_to = "mean_msic") |>
  mutate(species = recode(species, at = "A. thaliana", zs = "Z. stapffii"))


# fig: mean CDS MSIC against protein identity -----------------------------

p_pident <- ggplot(paired_long,
                   aes(x = pident, y = mean_msic, colour = species)) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.5) +
  geom_point(aes(shape = gene_set), size = 1.4) +
  facet_wrap(~model, nrow = 2) +
  scale_colour_manual(values = species_colors) +
  scale_shape_manual(values = c(BUSCO = 16, TPS = 17)) +
  labs(x = "Identity of best BLASTp segment (%)",
       y = "Mean MSIC over CDS per transcript",
       colour = "Species", shape = "Gene set") +
  guides(colour = guide_legend(
    theme = theme(legend.text = element_text(face = "italic"))
  )) +
  theme(legend.position = "inside", legend.position.inside = c(0.84, 0.22))

ggsave(file.path(fig_dir, "msic_x_blastp_cds.pdf"), p_pident,
       width = 16, height = 10, units = "cm", device = cairo_pdf)
