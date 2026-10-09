# Letters and Addresses
# A FASTA holds letters; a GFF holds addresses

library(Biostrings)

toy_seqs <- DNAStringSet(c(
  GeneA = "ATGGCTTCAGGT",
  GeneB = "ATGAAACCC"
))

toy_seqs
width(toy_seqs)
names(toy_seqs)

reverseComplement(toy_seqs)
subseq(toy_seqs, 1, 3)

writeXStringSet(toy_seqs, "toy.fa")
back_in <- readDNAStringSet("toy.fa")


# addresses ---------------------------------------------------------------

toy_ranges <- IRanges(start = c(1, 7),
                      end = c(3, 12),
                      names = c("exon1", "exon2"))

toy_ranges
start(toy_ranges)
end(toy_ranges)
width(toy_ranges) # 1-based: end - start + 1

GeneA <- toy_seqs[["GeneA"]]
extractAt(GeneA, toy_ranges)

intron_range <- gaps(toy_ranges)
extractAt(toy_seqs, at = intron_range)


# -------------------------------------------------------------------------

library(GenomicRanges)
library(Rsamtools)

indexFa("toy.fa")
toy_genome <- FaFile("toy.fa")

toy_gr <- GRanges(
  seqnames = c("GeneA", "GeneA", "GeneB"),
  ranges = IRanges(start = c(1, 7, 4), end = c(3, 12, 6)),
  strand = c("+", "+", "-")
)

toy_gr
getSeq(toy_genome, toy_gr)


# address container: GFF --------------------------------------------------

# gene -> mRNA -> two exons, all on sequence "geneA"
toy_gff_lines <- c(
  "##gff-version 3",
  "GeneA\ttoy\tgene\t1\t12\t.\t+\t.\tID=gene1",
  "GeneA\ttoy\tmRNA\t1\t12\t.\t+\t.\tID=tx1;Parent=gene1",
  "GeneA\ttoy\texon\t1\t3\t.\t+\t.\tID=ex1;Parent=tx1",
  "GeneA\ttoy\texon\t7\t12\t.\t+\t.\tID=ex2;Parent=tx1"
)
writeLines(toy_gff_lines, "toy.gff3")

library(rtracklayer)
toy_gff <- import("toy.gff3")

toy_gff

exons_geneA <- toy_gff[toy_gff$type == "exon"]

getSeq(toy_genome, exons_geneA)


# -------------------------------------------------------------------------

library(txdbmaker)
library(GenomicRanges)

toy_txdb <- makeTxDbFromGFF("toy.gff3", format = "gff3")

genes(toy_txdb)
exons_by_tx <- exonsBy(toy_txdb, by = "tx", use.names = TRUE)
introns_by_tx <- intronsByTranscript(toy_txdb, use.names = TRUE)

exons_by_tx
introns_by_tx



# -------------------------------------------------------------------------

library(rtracklayer)
library(Rsamtools)
library(Biostrings)

setwd("~/plantcad-zygophyllum-orthologs/")

genomes <- Rsamtools::FaFile("data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa")
gff <- rtracklayer::import("data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3")
genes <- gff[gff$type == "gene"]

# how many genes are on each chromosome (=seqnames)
table(seqnames(genes))

target <- genes[genes$gene_id == "AT1G71220"]
as.character(seqnames(target))
start(target)
end(target)
as.character(strand(target))
width(target)

# Which is the longest gene in the genome, and how long is it?
genes[which.max(width(genes))]$gene_id
which.max(width(genes))

# Cut out the sequence of AT1G71220 and compute its GC
target <- genes[genes$gene_id == "AT1G71220"]
target_seq <- getSeq(genomes, target)
letterFrequency(target_seq, letters = "GC", as.prob = T)

# Write a FASTA of all minus-strand genes on chromosome 2 longer
#  than 5,000 bp, with gene IDs as headers. This one uses all three steps.
targets <- genes[
  strand(genes) == "-" &
    seqnames(genes) == "2" &
    width(genes) > 5000
]

targets_seq <- getSeq(genomes, targets)
# headers
names(targets_seq) <- targets$gene_id
writeXStringSet(targets_seq, "learn_out.fa")



# GRanges list and split() ------------------------------------------------

library(GenomicRanges)

# A flat GRanges: six exons from three genes
exons <- GRanges(
  seqnames = c("chr1", "chr1", "chr1", "chr1", "chr1", "chr2"),
  ranges   = IRanges(start = c(100, 300, 600, 5000, 5800, 1000),
                     end   = c(200, 400, 700, 5300, 6000, 1200)),
  strand   = c("+", "+", "+", "-", "-", "+"),
  gene_id  = c("geneA", "geneA", "geneA", "geneB", "geneB", "geneC"),
  exon_id  = paste0("ex", 1:6)
)

# Group the exons by gene
exons_by_gene <- split(exons, exons$gene_id)
exons_by_gene

length(exons_by_gene)
length(exons_by_gene$geneA)
length(exons_by_gene$geneB)

# most funcs work on each element of list
width(exons_by_gene)
width(exons_by_gene$geneA)

sum(width(exons_by_gene))
start(exons_by_gene)
end(exons_by_gene)

range(exons_by_gene)


unlist(exons_by_gene)

# exercises
strand(exons_by_gene)

# What the exercise intended: flat GRanges, flat strand vector
exons_by_strand <- split(exons, strand(exons))
exons_by_strand



# -------------------------------------------------------------------------

toy <- GRanges(
  seqnames = c("chrA", "chrB"),
  ranges   = IRanges(start = c(10, 50), width = 20)
)
seqlevels(toy)
#> [1] "chrA" "chrB"

# Lengths for more contigs than toy has, in a different order
all_lengths <- c(chrB = 60, chrC = 500, chrA = 100)

seqlengths(toy) <- all_lengths[seqlevels(toy)]
seqinfo(toy)
#> chrA 100, chrB 60
