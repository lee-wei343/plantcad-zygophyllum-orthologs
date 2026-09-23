library(Rsamtools)
library(rtracklayer)

setwd("~/plantcad-zygophyllum-orthologs/")


# loading -----------------------------------------------------------------

genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
if (!file.exists(paste0(genome_path, ".fai"))) {
  indexFa(genome_path)
}
genome_fa <- FaFile(genome_path)


gff <- import("data/Z_stapfii/braker.gff3")
genes <- gff[gff$type == "gene"]

conserved <- data.table::fread("out/busco/conserved.tsv")

# quick look --------------------------------------------------------------

contig_lengths <- seqlengths(seqinfo(genome_fa))

# n contigs
length(contig_lengths)

which.max(contig_lengths)
max(contig_lengths)

table(gff$type)

exon_rows <- gff[gff$type == "exon"]
head(mcols(exon_rows))

# BUSCO cords to GRanges --------------------------------------------------

# GRanges always start as lower, BUSCO gave higher as start when minus strand
hit_start <- pmin(conserved$gene_start, conserved$gene_end)
hit_end <- pmax(conserved$gene_start, conserved$gene_end)

busco_hits <- GRanges(
  seqnames = conserved$contig,
  ranges = IRanges(start = hit_start, end = hit_end),
  strand = conserved$strand
)

# BUSCO cord to Braker gene -----------------------------------------------

# start and end of BUSCO and Braker might differ, but the will overlap
hits <- findOverlaps(busco_hits, genes)

# Check if sussy
# Each BUSCO matches one Braker
# It does not:
# matches no Braker;
# overlapps two braker genes;
# two BUSCO match to one Braker gene.
stopifnot(
  length(busco_hits) == length(hits),
  !anyDuplicated(queryHits(hits)),
  !anyDuplicated(subjectHits(hits)),
  identical(queryHits(hits), seq_along(busco_hits))
)

conserved$braker_gene_id <- genes$ID[subjectHits(hits)]

#conserved[, c("locus", "braker_gene_id")]
