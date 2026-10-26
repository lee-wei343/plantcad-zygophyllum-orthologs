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
blastn_dir <- file.path(out_dir, "blastn")
blastp_dir <- file.path(out_dir, "blastp")
plot_dir <- file.path(out_dir, "dot_plot")
dir.create(blastn_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(blastp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(plot_dir, showWarnings = FALSE, recursive = TRUE)

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


# make x_tx_by_cds for blastp ---------------------------------------------

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

# functions ---------------------------------------------------------------

process_pair <- function(agic_tx, bgc_tx) {
  pair_name <- paste0(agic_tx, "_with_", bgc_tx)
  message("processing ", pair_name)
  
  blast_dotplot(
    at_genome_fa = at_genome$fasta,
    zs_genome_fa = zs_genome$fasta,
    at_tx_list = at_tx_list,
    zs_tx_list = zs_tx_list,
    agic_tx = agic_tx,
    bgc_tx = bgc_tx,
    type = "blastn",
    outdir = blastn_dir,
    plot_dir = plot_dir
  )
  
  
  blast_dotplot(
    at_genome_fa = at_genome$fasta,
    zs_genome_fa = zs_genome$fasta,
    at_tx_list = at_tx_list,
    zs_tx_list = zs_tx_list,
    agic_tx = agic_tx,
    bgc_tx = bgc_tx,
    type = "blastp",
    outdir = blastp_dir,
    plot_dir = plot_dir,
    at_cds_by_tx = at_cds_by_tx,
    zs_cds_by_tx = zs_cds_by_tx
  )
  
}


# run ---------------------------------------------------------------------

walk2(tx_pair_keys$agic_tx, tx_pair_keys$bgc_tx, process_pair)
