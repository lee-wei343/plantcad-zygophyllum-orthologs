library(rtracklayer)
library(Rsamtools)
setwd("~/plantcad-zygophyllum-orthologs/")

get_exon_marked_seq <- function(gene_id, genome) {
  tx_id <- paste0(gene_id, ".t1")
  exons_current_gene <- exons_by_tx[[tx_id]]

  tx_span <- range(exons_current_gene)
  tx_span_flanked <- resize(
    tx_span,
    width = width(tx_span) + 1000,
    fix = "end"
  )
  tx_span_flanked <- resize(
    tx_span_flanked,
    width = width(tx_span_flanked) + 1000,
    fix = "start"
  )

  # contigs are far shorter than chromosomes, so keep the flank inside the contig
  seqlengths(tx_span_flanked) <- contig_lengths[seqlevels(tx_span_flanked)]
  tx_span_flanked <- trim(tx_span_flanked)

  if (as.character(strand(tx_span_flanked)) == "+") {
    current_strand <- "+"
    tx_span_start <- start(tx_span_flanked)
    exon_from <- start(exons_current_gene) - tx_span_start + 1
    exon_to <- end(exons_current_gene) - tx_span_start + 1
  } else {
    current_strand <- "-"
    tx_span_end <- end(tx_span_flanked)
    exon_to <- tx_span_end - start(exons_current_gene) + 1
    exon_from <- tx_span_end - end(exons_current_gene) + 1
  }

  chromosome <- as.character(seqnames(tx_span_flanked))

  flanked_seq <- getSeq(genome, tx_span_flanked)

  exonic_flanked_seq <- tolower(as.character(flanked_seq))

  # toupper EXONS
  for (i in seq_along(exons_current_gene)) {
    substr(exonic_flanked_seq, exon_from[i], exon_to[i]) <- toupper(substr(exonic_flanked_seq, exon_from[i], exon_to[i]))
  }

  header <- paste0(
    ">", tx_id, "|",
    gene_id, "|",
    chromosome, "|",
    current_strand, "|",
    start(tx_span_flanked), "-", end(tx_span_flanked), "|",
    "flank:1000"
  )

  return(list(
    header = header,
    exonic_flanked_seq = exonic_flanked_seq
  ))
}


# import conserved target genes
conserved <- data.table::fread("out/busco/conserved.tsv")

genome_path <- "data/Z_stapfii/Z_stapfii_genome.fa"
if (!file.exists(paste0(genome_path, ".fai"))) {
  indexFa(genome_path)
}
genome <- FaFile(genome_path)
contig_lengths <- seqlengths(seqinfo(genome))

gff_path <- "data/Z_stapfii/braker.gff3"
gff <- import(gff_path)

# braker rows carry transcript_id directly, but makeTxDbFromGFF() doesn't work :(
exons <- gff[gff$type == "exon"]
exons_by_tx <- split(exons, exons$transcript_id)

genes <- gff[gff$type == "gene"]

# conserved.tsv carries the A. thaliana locus plus the Z. stapfii BUSCO hit
# coordinates - there is no braker gene id in it, so resolve the gene by
# overlap. BUSCO reports CDS bounds and the gff reports gene bounds, so the
# two differ by a few bp and cannot be matched by equality.
busco_hits <- GRanges(
  seqnames = conserved$contig,
  ranges = IRanges(
    start = pmin(conserved$gene_start, conserved$gene_end),
    end = pmax(conserved$gene_start, conserved$gene_end)
  ),
  strand = conserved$strand
)

hits <- findOverlaps(busco_hits, genes)

# fail before writing anything rather than leaving a partial fasta behind:
# every BUSCO hit must resolve to exactly one gene and no two hits may share
# a gene, otherwise a fingerprint would carry the wrong A. thaliana locus
stopifnot(
  length(hits) == length(busco_hits),
  !anyDuplicated(queryHits(hits)),
  !anyDuplicated(subjectHits(hits))
)

conserved$gene_id <- as.character(genes$ID)[subjectHits(hits)]
conserved_genes <- conserved$gene_id

# braker annotates some genes with a second isoform; .t1 is taken throughout,
# mirroring the .1 convention used for A. thaliana
stopifnot(all(paste0(conserved_genes, ".t1") %in% names(exons_by_tx)))

if (file.exists("out/gene_region_fingerprint.Z_stapfii.fa")) {
  stop("file 'gene_region_fingerprint.Z_stapfii.fa' already exists!")
}
for (gene in conserved_genes) {
  res <- get_exon_marked_seq(gene_id = gene, genome = genome)
  cat(res$header,
    res$exonic_flanked_seq,
    sep = "\n",
    file = "out/gene_region_fingerprint.Z_stapfii.fa",
    append = TRUE
  )
}
