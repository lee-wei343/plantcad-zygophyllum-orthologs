library(readr)
library(rtracklayer)
library(Rsamtools)
library(txdbmaker)



setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")

at_gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"
zs_gff_path <- "data/Z_stapfii/braker.gff3"

flank_bp <- 1000L

# load --------------------------------------------------------------------

tps <- read_tsv("05-additional-genes/tps_orthos_in_busco.tsv")


Atha_genome <- load_genome(
  genome_path = "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
  gff_path = at_gff_path
)

Zsta_genome <- load_genome(
  genome_path = "data/Z_stapfii/Z_stapfii_genome.fa",
  gff_path = zs_gff_path
)

at_exons_by_tx <- get_exons_by_tx(
  gff_path = at_gff_path
)

at_exons_by_tx <- at_exons_by_tx[
  names(at_exons_by_tx) %in% paste0(tps$Atha, ".1")
]

zs_exons_by_tx <- split_exons_by_tx(
  gff = Zsta_genome$gff
)

zs_exons_by_tx <- zs_exons_by_tx[
  names(zs_exons_by_tx) %in% paste0(tps$Zsta, ".t1")
]

at_exons_by_tx_fl <- add_flanks(
  exons_by_tx = at_exons_by_tx,
  flank_bp = flank_bp,
  genome_fa = Atha_genome$fasta
)

zs_exons_by_tx_fl <- add_flanks(
  exons_by_tx = at_exons_by_tx,
  flank_bp = flank_bp,
  genome_fa = Atha_genome$fasta
)

