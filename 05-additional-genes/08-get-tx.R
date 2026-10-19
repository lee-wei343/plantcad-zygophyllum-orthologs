library(readr)
library(rtracklayer)
library(Rsamtools)


setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")


# load --------------------------------------------------------------------

tps <- read_tsv("05-additional-genes/tps_orthos_in_busco.tsv")

Zsta_genome <- load_genome(
  genome_path = "data/Z_stapfii/Z_stapfii_genome.fa",
  gff_path = "data/Z_stapfii/braker.gff3"
)
