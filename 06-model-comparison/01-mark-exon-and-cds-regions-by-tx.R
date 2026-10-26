library(tidyverse)
library(stringr)
library(Rsamtools)
library(txdbmaker)
library(rtracklayer)
library(GenomicRanges)
library(Biostrings)


# setup -------------------------------------------------------------------

setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")
source("04-top15/R/plot-helpers.R")

out_dir <- "06-model-comparison"
if (!dir.exists(out_dir)) {
  dir.create(out_dir, showWarnings = FALSE)
}

at_genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
at_gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"

zs_genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
zs_gff_path <- "data/Z_stapfii/braker.gff3"

flank_bp <- 1000L


# functions ---------------------------------------------------------------

mark_positions <- function(span, features) {
  if (is.null(features)) {
    return(rep(FALSE, width(span)))
  }
  rel <- shift(ranges(features), 1L - start(span))
  in_feature <- as.vector(coverage(rel, width = width(span)) > 0)
  if (as.character(strand(span)) == "-") {
    in_feature <- rev(in_feature)
  }
  in_feature
}

mark_positions_by_tx <- function(tx_pair_keys, at_genome, zs_genome) {
  
  # make x_tx_by_cds --------------------------------------------------------
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
  at_cds_by_tx <- at_cds_by_tx[names(at_cds_by_tx) %in% tx_pair_keys$agic_tx]
  zs_cds_by_tx <- cdsBy(zs_tx_db, by = "tx", use.names = TRUE)
  zs_cds_by_tx <- zs_cds_by_tx[names(zs_cds_by_tx) %in% tx_pair_keys$bgc_tx]
  
  
  # make x_tx_by_exons ------------------------------------------------------
  
  at_exons_by_tx <- exonsBy(at_tx_db, by = "tx", use.names = TRUE)
  at_exons_by_tx <- at_exons_by_tx[names(at_exons_by_tx) %in% tx_pair_keys$agic_tx]
  zs_exons_by_tx <- exonsBy(zs_tx_db, by = "tx", use.names = TRUE)
  zs_exons_by_tx <- zs_exons_by_tx[names(zs_exons_by_tx) %in% tx_pair_keys$bgc_tx]
  
  
  
  # get spans by tx ---------------------------------------------------------
  
  at_exon_spans <- add_flanks(at_exons_by_tx, flank_bp, at_genome$fasta)
  zs_exon_spans <- add_flanks(zs_exons_by_tx, flank_bp, zs_genome$fasta)
  
  
  
  
  at_pos_df <- map(
    set_names(names(at_exon_spans)),
    \(x) {
      span <- at_exon_spans[x]
      tibble(
        pos = seq_len(width(span)),
        in_exon = mark_positions(span, at_exons_by_tx[[x]]),
        in_cds = mark_positions(span, at_cds_by_tx[x][[1]])
      )
    }
  ) |>
    list_rbind(names_to = "tx_id")
  
  
  zs_pos_df <- map(
    set_names(names(zs_exon_spans)),
    \(x) {
      span <- zs_exon_spans[x]
      tibble(
        pos = seq_len(width(span)),
        in_exon = mark_positions(span, zs_exons_by_tx[[x]]),
        in_cds = mark_positions(span, zs_cds_by_tx[[x]])
      )
    }
  ) |>
    list_rbind(names_to = "tx_id")
  
  # at_pos_long <- bind_rows(at_pos_df, .id = "tx_id")
  # zs_pos_long <- bind_rows(zs_pos_df, .id = "tx_id")
   
  list(at_pos = at_pos_df, zs_pos = zs_pos_df)
} 

# load --------------------------------------------------------------------

tx_pair_keys_tps <- read_tsv(
  "05-additional-genes/tps_orthos_in_busco.tsv"
) |>
  mutate(
    agic_tx = paste0(Atha, ".1"),
    bgc_tx = paste0(Zsta, ".t1")
  ) |>
  dplyr::select(agic_tx, bgc_tx)

tx_pair_keys_top15 <- read_tsv("04-top15/out/busco_top15.tsv") |>
  mutate(agic_tx = paste0(agic, ".1"),
         bgc_tx = paste0(braker_gene_id, ".t1")) |>
  dplyr::select(agic_tx, bgc_tx)

at_genome <- load_genome(
  genome_path = at_genome_path,
  gff_path = at_gff_path
)

zs_genome <- load_genome(
  genome_path = zs_genome_path,
  gff_path = zs_gff_path
)

tps_marks <- mark_positions_by_tx(
  tx_pair_keys = tx_pair_keys_tps,
  at_genome = at_genome,
  zs_genome = zs_genome
)

top15_marks <- mark_positions_by_tx(
  tx_pair_keys = tx_pair_keys_top15,
  at_genome = at_genome,
  zs_genome = zs_genome
) 
  
saveRDS(tps_marks, file = "06-model-comparison/tps_marks.rds")
saveRDS(top15_marks, file = "06-model-comparison/top15_marks.rds")
