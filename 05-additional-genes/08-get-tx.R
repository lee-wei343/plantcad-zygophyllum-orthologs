library(readr)
library(rtracklayer)
library(Rsamtools)


setwd("~/plantcad-zygophyllum-orthologs/")


# load --------------------------------------------------------------------

tps <- read_tsv("05-additional-genes/tps_orthos_in_busco.tsv")

Zsta_genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
if (!file.exists(paste0(Zsta_genome_path, ".fai"))) {
  indexFa(Zsta_genome_path)
}
Zsta_genome_fa <- FaFile(Zsta_genome_path)
