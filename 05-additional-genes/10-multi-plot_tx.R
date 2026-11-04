library(zoo)
library(tidyverse)
library(stringr)
library(patchwork)
library(Rsamtools)
library(txdbmaker)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)

setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")
source("04-top15/R/plot-helpers.R")


# settings ----------------------------------------------------------------

out_dir <- "05-additional-genes"
multi_path <- file.path(out_dir, "out", "multi_plot")
dir.create(multi_path, showWarnings = FALSE, recursive = TRUE)

at_genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
at_gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"

zs_genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
zs_gff_path <- "data/Z_stapfii/braker.gff3"

flank_bp <- 1000L


# color blind friendly palette --------------------------------------------

palette <- palette.colors(palette = "Okabe-Ito", names = TRUE)
class_colors <- c(
  "low" = palette[["blue"]],
  "neutral" = palette[["gray"]],
  "high" = palette[["vermillion"]]
)


# load --------------------------------------------------------------------

tx_pair_keys <- read_tsv(
  "05-additional-genes/tps_orthos_in_busco.tsv"
  ) |>
  mutate(
    agic_tx = paste0(Atha, ".1"),
    bgc_tx = paste0(Zsta, ".t1")
  ) |>
  dplyr::select(agic_tx, bgc_tx)


at_genome <- load_genome(
  genome_path = at_genome_path,
  gff_path = at_gff_path
)

zs_genome <- load_genome(
  genome_path = zs_genome_path,
  gff_path = zs_gff_path
)


# make x_tx_by_cds for  --------------------------------------------------

at_tx_db <- makeTxDbFromGFF(
  at_gff_path,
  format = "gff3",
  )

# TxDb of the BRAKER gff: transcript names come from Name, not transcript_id
zs_gff <- zs_genome$gff
zs_gff$transcript_id <- NULL
zs_gff$Name <- ifelse(zs_gff$type == "mRNA", zs_gff$ID, NA_character_)
zs_tx_db <- makeTxDbFromGRanges(zs_gff)

at_cds_by_tx <- cdsBy(at_tx_db, by = "tx", use.names = TRUE)
zs_cds_by_tx <- cdsBy(zs_tx_db, by = "tx", use.names = TRUE)


# tx lists ----------------------------------------------------------------

at_exons_by_tx <- exonsBy(at_tx_db, by = "tx", use.names = TRUE)
at_exons_by_tx <- at_exons_by_tx[names(at_exons_by_tx) %in% tx_pair_keys$agic_tx]
at_exon_spans <- add_flanks(at_exons_by_tx, flank_bp, at_genome$fasta)

at_tx_list <- make_tx_list(
  exon_spans = at_exon_spans,
  exons_by_tx = at_exons_by_tx,
  headers = make_header(at_exon_spans, at_exons_by_tx, flank_bp)
)

zs_exons_by_tx <- split_exons_by_tx(zs_genome$gff)
zs_exons_by_tx <- zs_exons_by_tx[names(zs_exons_by_tx) %in% tx_pair_keys$bgc_tx]
zs_exon_spans <- add_flanks(zs_exons_by_tx, flank_bp, zs_genome$fasta)

zs_tx_list <- make_tx_list(
  exon_spans = zs_exon_spans,
  exons_by_tx = zs_exons_by_tx,
  headers = make_header(zs_exon_spans, zs_exons_by_tx, flank_bp, seq_label = "ctg")
)

stopifnot(
  all(tx_pair_keys$agic_tx %in% names(at_tx_list)),
  all(tx_pair_keys$bgc_tx %in% names(zs_tx_list))
)


# load --------------------------------------------------------------------

results_dir <- function(species, model) {
  file.path(out_dir, "out", paste0(species, "_tps_", model, "_results"))
}

# PlantCAD1 / PlantCAD2-small: LLMPipeGUI, one csv per transcript
read_msic_llmpipe <- function(species, model, tx_id) {
  dir <- results_dir(species, model)
  read_msic(
    results_dir = dir,
    file_prefix = paste0(species, "_tps_", model, "_scores_table"),
    index = make_tx_index(readLines(file.path(dir, "fileIDs.txt"))),
    tx_id = tx_id
  )
}

# PlantCAD2-large / Botanic1: one msic.csv for all transcripts
read_msic_combined <- function(species, model, tx_id) {
  read_csv(
    file.path(results_dir(species, model), "msic.csv"),
    show_col_types = FALSE
  ) |>
    mutate(tx = str_extract(header, "^[^|]+")) |>
    filter(tx == tx_id) |>
    dplyr::rename(MSIC = msic) |>
    add_pos() |>
    add_class() |>
    add_rolling_mean()
}

read_all_models <- function(species, tx_id) {
  list(
    PCAD1 = read_msic_llmpipe(species, "PCAD1", tx_id),
    PCAD2S = read_msic_llmpipe(species, "PCAD2S", tx_id),
    PCAD2L = read_msic_combined(species, "PCAD2L", tx_id),
    BOT1S = read_msic_combined(species, "Botanic1S", tx_id),
    BOT1L = read_msic_combined(species, "Botanic1L", tx_id)
  )
}

# every model must have scored the same bases, or the facets are misaligned
check_same_bases <- function(msic_table, expected_length) {
  refs <- map(msic_table, \(df) tolower(df$ref))
  stopifnot(
    all(lengths(refs) == expected_length),
    all(map_lgl(refs, \(ref) identical(ref, refs[[1]])))
  )
}


# -------------------------------------------------------------------------

for (i in seq_len(nrow(tx_pair_keys))) {

  current_agic_tx <- tx_pair_keys$agic_tx[i]
  current_bgc_tx <- tx_pair_keys$bgc_tx[i]

  msic_table_at <- read_all_models("a_thaliana", current_agic_tx)
  msic_table_zs <- read_all_models("z_stapffii", current_bgc_tx)

  check_same_bases(msic_table_at, at_tx_list[[current_agic_tx]]$end)
  check_same_bases(msic_table_zs, zs_tx_list[[current_bgc_tx]]$end)
  

  msic_long_at <- bind_rows(msic_table_at, .id = "run") |>
    mutate(run = factor(run, levels = names(msic_table_at)))

  msic_long_zs <- bind_rows(msic_table_zs, .id = "run") |>
    mutate(run = factor(run, levels = names(msic_table_zs)))

  x_range  <- c(1, max(msic_long_at$pos, msic_long_zs$pos))
  x_scale  <- scale_x_continuous(expand = expansion(mult = 0.01))

  msic_plot_at <- ggplot(msic_long_at, aes(x = pos, y = MSIC)) +
    geom_col(aes(fill = class), width = 1L) +
    geom_hline(yintercept = 0) +
    geom_line(aes(y = rolling_mean), linewidth = 0.3) +
    scale_fill_manual(values = class_colors) +
    facet_grid(run ~ .) +
    x_scale +
    coord_cartesian(xlim = x_range, ylim = c(-1, 1)) +
    labs(y = "MSIC", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank())


  msic_plot_zs <- ggplot(msic_long_zs, aes(x = pos, y = MSIC)) +
    geom_col(aes(fill = class), width = 1L) +
    geom_hline(yintercept = 0) +
    geom_line(aes(y = rolling_mean), linewidth = 0.3) +
    scale_fill_manual(values = class_colors) +
    facet_grid(run ~ .) +
    x_scale +
    coord_cartesian(xlim = x_range, ylim = c(-1, 1)) +
    labs(y = "MSIC", x = NULL) +
    theme_minimal() +
    theme(axis.text.x = element_blank())

  current_at_tx_list <- at_tx_list[[current_agic_tx]]
  tx_plot_at <- plot_transcript_simple(
    start = current_at_tx_list$start,
    end = current_at_tx_list$end,
    strand = current_at_tx_list$strand,
    exon_starts = current_at_tx_list$relative_exon_start,
    exon_ends = current_at_tx_list$relative_exon_end,
    title = current_at_tx_list$header
  ) +
    x_scale +
    coord_cartesian(xlim = x_range)

  current_zs_tx_list <- zs_tx_list[[current_bgc_tx]]
  tx_plot_zs <- plot_transcript_simple(
    start = current_zs_tx_list$start,
    end = current_zs_tx_list$end,
    strand = current_zs_tx_list$strand,
    exon_starts = current_zs_tx_list$relative_exon_start,
    exon_ends = current_zs_tx_list$relative_exon_end,
    title = current_zs_tx_list$header
  ) +
    x_scale +
    coord_cartesian(xlim = x_range)

  combined_plot <- tx_plot_at / msic_plot_at / tx_plot_zs / msic_plot_zs +
    plot_layout(heights = c(1, 8, 1, 8), guides = "collect")

  combined_plot

  ggsave(paste0(multi_path, "/", current_agic_tx, "_with_", current_bgc_tx, ".png"),
         plot = combined_plot,
         dpi = 300,
         height = 7,
         width = 12)
}