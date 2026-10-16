library(stringr)
library(tibble)
library(purrr)
library(readr)
library(ape)

setwd("~/plantcad-zygophyllum-orthologs/")


# functions ---------------------------------------------------------------

make_tree <- function(lines, orthogroup) {
  og_line <- lines[startsWith(lines, orthogroup)]
  tree <- read.tree(text = str_replace(og_line, "^OG\\d+:\\s*", ""))
}

make_pruned_tree <- function(lines, orthogroup, prune) {
  tree <- make_tree(lines = lines, orthogroup = orthogroup)
  tip <- grep(pattern = prune, x = tree$tip.label, ignore.case = TRUE)
  pruned_tree <- drop.tip(tree, tip)
}

write_tree_plot <- function(tree, name, out_dir) {
  png(
    file.path(out_dir, paste0(name, "_tree.png")),
    width = 8,
    height = 5,
    units = "in",
    res = 300
    )
  plot(tree)
  dev.off()
}

# example -----------------------------------------------------------------

tr <- read.tree(
  text = "((Zsta_g0500,(TmoA08G003000,TmoB08G002950)),(Zsta_g0501,TmoB08G002960));"
)
plot(tr)



# -------------------------------------------------------------------------

out_dir <- file.path("05-additional-genes", "ogs_tree")
gene_tree_path <- file.path(
  "02-orthofinder",
  "orthofinder_longest_isoform",
  "Results_Oct15",
  "Resolved_Gene_Trees",
  "Resolved_Gene_Trees.txt"
  )

dir.create(out_dir, showWarnings = FALSE)


ogs <- readLines("05-additional-genes/ogs.txt")
lines <- readLines(gene_tree_path)

walk(ogs, \(og) {
  tree <- make_tree(lines = lines, orthogroup = og)
  write_tree_plot(tree = tree, name = og, out_dir = out_dir)
})

walk(ogs, \(og) {
  tree <- make_pruned_tree(
    lines = lines,
    orthogroup = og,
    prune = "^(o.*sativa*|a.*thaliana).*"
    )
  write_tree_plot(tree = tree, name = paste0(og, "_pruned"), out_dir = out_dir)
})
