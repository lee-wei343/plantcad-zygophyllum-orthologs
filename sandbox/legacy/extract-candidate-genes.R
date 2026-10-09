library(GenomicFeatures)
library(txdbmaker)

txdb <- makeTxDbFromGFF("data/Z_stapfii/braker.gff3")

gene_ranges <- genes(txdb)
gene_lengths <- width(gene_ranges)

exons_per_gene <- exonsBy(txdb, by = "gene")
exon_counts <- lengths(exons_per_gene)

length_ok <- names(gene_ranges)[gene_lengths >= 5000 & gene_lengths <= 10000]
exon_ok <- names(exon_counts)[exon_counts > 5]

selected_genes <- intersect(length_ok, exon_ok)

result <- gene_ranges[selected_genes]

dir.create("out", recursive = TRUE, showWarnings = FALSE)
writeLines(names(result), "out/selected_candidate-gene_ids.txt")
