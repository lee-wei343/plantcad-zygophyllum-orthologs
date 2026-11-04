library(tidyverse)
library(ggplot2)

setwd("~/plantcad-zygophyllum-orthologs/06-model-comparison/")

theme_set(theme_minimal())

df <- read_csv("msic_long.csv", show_col_types = FALSE)

thr <- 0.5

plot_dir <- "out"
dir.create(plot_dir, showWarnings = FALSE)


# PlantCAD1 ---------------------------------------------------------------

p1e <- df |> 
  filter(model == "PCAD1", in_exon == TRUE)

p1e_sup <- p1e |>
  filter(!is.na(msic)) |>
  mutate(region = cut(
    msic, c(-Inf, -thr, thr, Inf),
    labels = c("low", "mid", "high")
  )) |>
  count(species, is_busco, region) |>
  mutate(
    proportion = n / sum(n),
    label = scales::percent(proportion, accuracy = 0.1),
    x = c(low = -thr, mid = 0, high = thr)[as.character(region)],
    hjust = c(low = 1.1, mid = 0.5, high = -0.1)[as.character(region)],
    .by = c(species, is_busco)
  )

p1e_plot <- ggplot(data = p1e) +
  geom_histogram(
    aes(x = msic, y = after_stat(count / ave(count, PANEL, FUN = sum)))
    ) +
  geom_vline(xintercept = c(-thr, thr)) +
  geom_text(data = p1e_sup,
            aes(x = x, y = Inf, label = label, hjust = hjust),
            vjust = 1.3, size = 3) +
  facet_grid(species ~ is_busco,
             labeller = labeller(
               is_busco = c(`TRUE` = "BUSCO", `FALSE` = "TPS"),
               species = c("a_thaliana" = "A. thaliana",
                           "z_stapffii" = "Z. stapffii"),
               name = NULL),
             scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0.2))
    ) +
  labs(y = "relative frequency", x = "MSIC") +
  labs(title = "PCAD1 exons")

p1e_plot

ggsave(file.path(plot_dir, "msic_freq_PCAD1_exons.png"),
       plot = p1e_plot,
       dpi = 300,
       height = 5,
       width = 8)



# PlantCAD2S --------------------------------------------------------------

p2se <- df |> 
  filter(model == "PCAD2S", in_exon == TRUE)

p2se_sup <- p2se |>
  filter(!is.na(msic)) |>
  mutate(region = cut(
    msic, c(-Inf, -thr, thr, Inf),
    labels = c("low", "mid", "high")
  )) |>
  count(species, is_busco, region) |>
  mutate(
    proportion = n / sum(n),
    label = scales::percent(proportion, accuracy = 0.1),
    x = c(low = -thr, mid = 0, high = thr)[as.character(region)],
    hjust = c(low = 1.1, mid = 0.5, high = -0.1)[as.character(region)],
    .by = c(species, is_busco)
  )

p2se_plot <- ggplot(data = p2se) +
  geom_histogram(
    aes(x = msic, y = after_stat(count / ave(count, PANEL, FUN = sum)))
  ) +
  geom_vline(xintercept = c(-thr, thr)) +
  geom_text(data = p2se_sup,
            aes(x = x, y = Inf, label = label, hjust = hjust),
            vjust = 1.3, size = 3) +
  facet_grid(species ~ is_busco,
             labeller = labeller(
               is_busco = c(`TRUE` = "BUSCO", `FALSE` = "TPS"),
               species = c("a_thaliana" = "A. thaliana",
                           "z_stapffii" = "Z. stapffii"),
               name = NULL),
             scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(y = "relative frequency", x = "MSIC") +
  labs(title = "PCAD2S exons")

p2se_plot

ggsave(file.path(plot_dir, "msic_freq_PCAD2S_exons.png"),
       plot = p2se_plot,
       dpi = 300,
       height = 5,
       width = 8)




# PlantCAD2L --------------------------------------------------------------

p2le <- df |> 
  filter(model == "PCAD2L", in_exon == TRUE)

p2le_sup <- p2le |>
  filter(!is.na(msic)) |>
  mutate(region = cut(
    msic, c(-Inf, -thr, thr, Inf),
    labels = c("low", "mid", "high")
  )) |>
  count(species, is_busco, region) |>
  mutate(
    proportion = n / sum(n),
    label = scales::percent(proportion, accuracy = 0.1),
    x = c(low = -thr, mid = 0, high = thr)[as.character(region)],
    hjust = c(low = 1.1, mid = 0.5, high = -0.1)[as.character(region)],
    .by = c(species, is_busco)
  )

p2le_plot <- ggplot(data = p2le) +
  geom_histogram(
    aes(x = msic, y = after_stat(count / ave(count, PANEL, FUN = sum)))
  ) +
  geom_vline(xintercept = c(-thr, thr)) +
  geom_text(data = p2le_sup,
            aes(x = x, y = Inf, label = label, hjust = hjust),
            vjust = 1.3, size = 3) +
  facet_grid(species ~ is_busco,
             labeller = labeller(
               is_busco = c(`TRUE` = "BUSCO", `FALSE` = "TPS"),
               species = c("a_thaliana" = "A. thaliana",
                           "z_stapffii" = "Z. stapffii"),
               name = NULL),
             scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(y = "relative frequency", x = "MSIC") +
  labs(title = "PCAD2L exons")

p2le_plot

ggsave(file.path(plot_dir, "msic_freq_PCAD2L_exons.png"),
       plot = p2le_plot,
       dpi = 300,
       height = 5,
       width = 8)



# Botanic1S ---------------------------------------------------------------

b1se <- df |> 
  filter(model == "Botanic1S", in_exon == TRUE)

b1se_sup <- b1se |>
  filter(!is.na(msic)) |>
  mutate(region = cut(
    msic, c(-Inf, -thr, thr, Inf),
    labels = c("low", "mid", "high")
  )) |>
  count(species, is_busco, region) |>
  mutate(
    proportion = n / sum(n),
    label = scales::percent(proportion, accuracy = 0.1),
    x = c(low = -thr, mid = 0, high = thr)[as.character(region)],
    hjust = c(low = 1.1, mid = 0.5, high = -0.1)[as.character(region)],
    .by = c(species, is_busco)
  )

b1se_plot <- ggplot(data = b1se) +
  geom_histogram(
    aes(x = msic, y = after_stat(count / ave(count, PANEL, FUN = sum)))
  ) +
  geom_vline(xintercept = c(-thr, thr)) +
  geom_text(data = b1se_sup,
            aes(x = x, y = Inf, label = label, hjust = hjust),
            vjust = 1.3, size = 3) +
  facet_grid(species ~ is_busco,
             labeller = labeller(
               is_busco = c(`TRUE` = "BUSCO", `FALSE` = "TPS"),
               species = c("a_thaliana" = "A. thaliana",
                           "z_stapffii" = "Z. stapffii"),
               name = NULL),
             scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(y = "relative frequency", x = "MSIC") +
  labs(title = "Botanic1S exons")

b1se_plot

ggsave(file.path(plot_dir, "msic_freq_Botanic1S_exons.png"),
       plot = b1se_plot,
       dpi = 300,
       height = 5,
       width = 8)



# Botanic1L ---------------------------------------------------------------

b1le <- df |> 
  filter(model == "Botanic1L", in_exon == TRUE)

b1le_sup <- b1le |>
  filter(!is.na(msic)) |>
  mutate(region = cut(
    msic, c(-Inf, -thr, thr, Inf),
    labels = c("low", "mid", "high")
  )) |>
  count(species, is_busco, region) |>
  mutate(
    proportion = n / sum(n),
    label = scales::percent(proportion, accuracy = 0.1),
    x = c(low = -thr, mid = 0, high = thr)[as.character(region)],
    hjust = c(low = 1.1, mid = 0.5, high = -0.1)[as.character(region)],
    .by = c(species, is_busco)
  )

b1le_plot <- ggplot(data = b1le) +
  geom_histogram(
    aes(x = msic, y = after_stat(count / ave(count, PANEL, FUN = sum)))
  ) +
  geom_vline(xintercept = c(-thr, thr)) +
  geom_text(data = b1le_sup,
            aes(x = x, y = Inf, label = label, hjust = hjust),
            vjust = 1.3, size = 3) +
  facet_grid(species ~ is_busco,
             labeller = labeller(
               is_busco = c(`TRUE` = "BUSCO", `FALSE` = "TPS"),
               species = c("a_thaliana" = "A. thaliana",
                           "z_stapffii" = "Z. stapffii"),
               name = NULL),
             scales = "free_y") +
  scale_y_continuous(
    labels = scales::label_percent(),
    expand = expansion(mult = c(0, 0.2))
  ) +
  labs(y = "relative frequency", x = "MSIC") +
  labs(title = "Botanic1L exons")

b1le_plot

ggsave(file.path(plot_dir, "msic_freq_Botanic1L_exons.png"),
       plot = b1le_plot,
       dpi = 300,
       height = 5,
       width = 8)

