library(rtracklayer)
library(Rsamtools)
library(txdbmaker)
setwd("~/plantcad-zygophyllum-orthologs/")

get_exon_marked_seq <- function(gene_id, genome) {
  tx_id <- paste0(gene_id, ".1")
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
    "flank:1000")
  
  return(list(
    header = header,
    exonic_flanked_seq = exonic_flanked_seq
  ))
}


# import conserved target genes
conserved <- data.table::fread("out/busco/conserved.tsv")
wanted_ids <- conserved$locus

genome_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
genome <- FaFile(genome_path)

gff_path <- "data/A_thaliana/Arabidopsis_thaliana.TAIR10.63.gff3"
gff <- import(gff_path)

# create db, for exon to transcript parent link
txdb <- makeTxDbFromGFF(gff_path, format = "gff3")
exons_by_tx <- exonsBy(txdb, by = "tx", use.names = TRUE)

genes <- gff[gff$type == "gene"]
conserved_genes <- genes[genes$gene_id %in% wanted_ids]

if (file.exists("out/gene_region_fingerprint.fa")) {
  stop("file 'gene_region_fingerprint.fa' already exists!")
}
for (gene in conserved_genes$gene_id) {

  res <- get_exon_marked_seq(gene_id = gene, genome = genome)
  cat(res$header,
    res$exonic_flanked_seq,
    sep = "\n",
    file = "out/gene_region_fingerprint.fa",
    append = TRUE
  )
}
