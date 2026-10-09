library(rtracklayer)
library(Rsamtools)
library(Biostrings)

setwd("~/plantcad-zygophyllum-orthologs/")


# extract gene objects from gff -------------------------------------------

gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"

gff <- import(gff_path)
length(gff)
table(gff$type)

genes <- gff[gff$type == "gene"]
length(genes)

conserved <- data.table::fread("out/busco/conserved.tsv")
wanted_ids <- conserved$locus

conserved_genes <- genes[genes$gene_id %in% wanted_ids]
length(conserved_genes)
conserved_genes


# extract gff genes from genome -------------------------------------------

genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
genome <- FaFile(genome_path)

gene_seqs <- getSeq(genome, conserved_genes)

names(gene_seqs) <- conserved_genes$gene_id

writeXStringSet(gene_seqs, "out/busco/conserved_A_thaliana.fa")