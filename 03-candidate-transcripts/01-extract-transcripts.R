library(Rsamtools)
library(rtracklayer)
library(txdbmaker)
library(ggplot2)
library(pwalign)

setwd("~/plantcad-zygophyllum-orthologs/")
source("R/extract-helpers.R")


# setup -------------------------------------------------------------------

flank_bp <- 1000L
save_fasta <- TRUE
make_plots <- TRUE

at_genome_path <- file.path(
  "data",
  "A_thaliana",
  "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
  )
at_gff_path <- file.path(
  "data",
  "A_thaliana",
  "Arabidopsis_thaliana.TAIR10.63.gff3"
)

zs_genome_path <- file.path(
  "data",
  "Z_stapfii",
  "Z_stapfii_genome.fa"
)
zs_gff_path <- file.path(
  "data",
  "Z_stapfii",
  "braker.gff3"
)

busco <- readr::read_tsv("01-busco/out/conserved_go_annotated.tsv")

# Arabidopsis thaliana ----------------------------------------------------

# at = Arabidopsis thaliana
# Arabidopsis Genome Initiative Code


# load AT -----------------------------------------------------------------

at_genome <- load_genome(genome_path = at_genome_path)

at_exons_by_tx <- get_exons_by_tx(gff_path = at_gff_path)


# Filter AT ---------------------------------------------------------------

at_exons_by_busco_tx <- at_exons_by_tx[
  names(at_exons_by_tx) %in% paste0(busco$agic, ".1")
]


# get Seq with header -----------------------------------------------------

# GRanges
at_ex_by_tx_flank <- add_flanks(
  exons_by_tx = at_exons_by_busco_tx,
  flank_bp = flank_bp,
  genome_fa = at_genome$fasta
  )

# DNA-string set
at_tx_seqs <- getSeq(at_genome$fasta, at_ex_by_tx_flank)
stopifnot(
  length(at_ex_by_tx_flank) == length(at_tx_seqs),
  identical(width(at_ex_by_tx_flank), width(at_tx_seqs))
)

headers <- make_header(
  exon_spans = at_ex_by_tx_flank,
  exons_by_tx = at_exons_by_busco_tx,
  flank_bp = flank_bp
)
names(at_tx_seqs) <- headers


# write FASTA -------------------------------------------------------------

if (save_fasta) {
  writeXStringSet(
    at_tx_seqs,
    filepath = "03-candidate-transcripts/out/tx_a_thaliana.fa"
    )
}


# Plot Seqs ---------------------------------------------------------------

if (make_plots) {
  plot_tx_list(
    tx_list = make_tx_list(
      exon_spans = at_ex_by_tx_flank,
      exons_by_tx = at_exons_by_busco_tx,
      headers = headers
    ),
    out_dir = "03-candidate-transcripts/out/a_thaliana_tx_plots"
  )
}


# Z stapffii --------------------------------------------------------------

# zs = Z stapffii

zs_genome <- load_genome(
  genome_path = zs_genome_path,
  gff_path = zs_gff_path
  )

zs_genes <- zs_genome$gff[zs_genome$gff$type == "gene"]


# quick look --------------------------------------------------------------

# contig_lengths <- seqlengths(seqinfo(genome_fa))
#
# # n contigs
# length(contig_lengths)
#
# which.max(contig_lengths)
# max(contig_lengths)
#
# table(gff$type)


# BUSCO cords to GRanges --------------------------------------------------

# GRanges always start as lower, BUSCO gave higher as start when minus strand
hit_start <- pmin(busco$gene_start, busco$gene_end)
hit_end <- pmax(busco$gene_start, busco$gene_end)

busco_hits <- GRanges(
  seqnames = busco$contig,
  ranges = IRanges(start = hit_start, end = hit_end),
  strand = busco$strand
)

# BUSCO cord to Braker gene -----------------------------------------------

# start and end of BUSCO and Braker might differ, but they will overlap
hits <- findOverlaps(busco_hits, zs_genes)

# Check if suspicious
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

busco$braker_gene_id <- zs_genes$ID[subjectHits(hits)]

# get exons for each transcript -------------------------------------------

zs_exons_by_tx <- split_exons_by_tx(zs_genome$gff)

# I take only the first transcripts (.t1 not .t2 etc.) of a gene
zs_tx_ids <- paste0(busco$braker_gene_id, ".t1")
stopifnot(all(zs_tx_ids %in% names(zs_exons_by_tx)))

zs_exons_by_busco_tx <- zs_exons_by_tx[zs_tx_ids]


# get Seq with header -----------------------------------------------------

# I could also allow trimming, and note it in FASTA header
zs_ex_by_tx_flank <- add_flanks(
  exons_by_tx = zs_exons_by_busco_tx,
  flank_bp = flank_bp,
  genome_fa = zs_genome$fasta
)
stopifnot(identical(names(zs_ex_by_tx_flank), zs_tx_ids))

zs_tx_seqs <- getSeq(zs_genome$fasta, zs_ex_by_tx_flank)
stopifnot(
  length(zs_tx_seqs) == length(zs_ex_by_tx_flank),
  identical(width(zs_tx_seqs), width(zs_ex_by_tx_flank))
)

# zs_tx_ids is in busco row order, so the AGI loci line up
zs_headers <- make_header(
  exon_spans = zs_ex_by_tx_flank,
  exons_by_tx = zs_exons_by_busco_tx,
  flank_bp = flank_bp,
  seq_label = "ctg",
  extra = busco$agic
)
names(zs_tx_seqs) <- zs_headers

if (save_fasta) {
  writeXStringSet(zs_tx_seqs, "03-candidate-transcripts/out/tx_z_stapffii.fa")
}

if (make_plots) {
  plot_tx_list(
    tx_list = make_tx_list(
      exon_spans = zs_ex_by_tx_flank,
      exons_by_tx = zs_exons_by_busco_tx,
      headers = zs_headers
    ),
    out_dir = "03-candidate-transcripts/out/z_stapfii_tx_plots"
  )
}


# Pairwise align ----------------------------------------------------------

# Pairwise global alignment
dna_matrix <- nucleotideSubstitutionMatrix(
  match = 2,
  mismatch = -1,
  baseOnly = TRUE
)

pw_align_scores <- numeric(nrow(busco))
pid_align_scores <- numeric(nrow(busco))
for (i in seq_len(nrow(busco))) {
  at_tx <- at_ex_by_tx_flank[paste0(busco$agic[i], ".1")]
  at_tx_seq <- getSeq(at_genome$fasta, at_tx)

  zs_tx <- zs_ex_by_tx_flank[paste0(busco$braker_gene_id[i], ".t1")]
  zs_tx_seq <- getSeq(zs_genome$fasta, zs_tx)

  global_align <- pairwiseAlignment(
    at_tx_seq,
    zs_tx_seq,
    type = "global",
    substitutionMatrix = dna_matrix,
    gapOpening = 5, gapExtension = 1
  )
  pw_align_scores[i] <- score(global_align)
  pid_align_scores[i] <- pid(global_align)
}

busco$pw_align_scores <- pw_align_scores
busco$pid_align_scores <- pid_align_scores


readr::write_tsv(busco, file = "03-candidate-transcripts/out/conserved_go_annotated_braker.tsv")
