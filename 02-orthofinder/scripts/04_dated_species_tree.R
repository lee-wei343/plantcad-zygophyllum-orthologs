library(ape)

setwd("~/plantcad-zygophyllum-orthologs/")

# Divergence times in Mya. Branch lengths sum to these node ages.
#   TimeTree pairwise API, precomputed_age (CI), queried 2026-11-04:
#     A. thaliana - P. trichocarpa   108 (104-111)
#     A. thaliana - V. vinifera      118 (115-121)
#     A. thaliana - O. sativa        160 (154-166)
#     A. thaliana - P. patens        483 (465-500)
#   van der Merwe et al. 2026, bioRxiv (preprint), BEAST, Fig. 6:
#     Tetraena - Z. stapffii group   20.76 (12.06-29.01)
#   Liu et al. 2023, Fig. 2C agrees with the 108 node:
#     T. mongolica - other rosids    100-108

out_dir <- "02-orthofinder/out"


# functions ---------------------------------------------------------------

plot_dated_tree <- function(tree, node_labels, edge_labels, file) {
  tree$tip.label <- gsub("_", " ", tree$tip.label)

  png(file, width = 11, height = 6, units = "in", res = 300)
  on.exit(dev.off())

  par(mar = c(4.5, 1, 1, 1))
  plot(tree, font = 3, cex = 1.05, label.offset = 6, no.margin = FALSE,
       x.lim = c(0, 640), edge.width = 1.6)

  for (label in node_labels) {
    node <- if (length(label$tips) == 2) {
      getMRCA(tree, label$tips)
    } else {
      Ntip(tree) + 1
    }
    nodelabels(label$text, node, frame = "none", cex = 0.9,
               adj = label$adj, font = 2)
  }

  for (label in edge_labels) {
    child <- if (length(label$tips) == 1) {
      which(tree$tip.label == label$tips)
    } else {
      getMRCA(tree, label$tips)
    }
    edgelabels(label$text, which(tree$edge[, 2] == child), frame = "none",
               cex = 0.85, adj = c(0.5, -0.5), col = "grey30")
  }

  axisPhylo(backward = TRUE)
  mtext("Million years ago", side = 1, line = 2.8)
}


# trees -------------------------------------------------------------------

# exact topology of the whiteboard: rosids is one polytomy
whiteboard <- read.tree(text = paste0(
  "(Physcomitrella_patens:483,",
  "(Oryza_sativa:160,",
  "(Vitis_vinifera:118,Populus_trichocarpa:118,",
  "(Tetraena_mongolica:20.8,Zygophyllum_stapffii:20.8):97.2,",
  "Arabidopsis_thaliana:118):42):323);"
))

# the species used in the thesis, without V. vinifera and P. trichocarpa
dated <- read.tree(text = paste0(
  "(Physcomitrella_patens:483,",
  "(Oryza_sativa:160,",
  "((Tetraena_mongolica:20.8,Zygophyllum_stapffii:20.8):87.2,",
  "Arabidopsis_thaliana:108):52):323);"
))


# plot --------------------------------------------------------------------

shared_nodes <- list(
  list(text = "Embryophytes  483", tips = "root", adj = c(-0.04, -0.6)),
  list(text = "angiosperms  160",
       tips = c("Oryza sativa", "Arabidopsis thaliana"), adj = c(1.04, -0.6)),
  list(text = "20.8",
       tips = c("Tetraena mongolica", "Zygophyllum stapffii"),
       adj = c(-0.2, 0.5))
)

whiteboard_nodes <- c(shared_nodes, list(
  list(text = "rosids  118",
       tips = c("Vitis vinifera", "Arabidopsis thaliana"), adj = c(1.04, -0.6))
))

dated_nodes <- c(shared_nodes, list(
  list(text = "rosids  108",
       tips = c("Zygophyllum stapffii", "Arabidopsis thaliana"),
       adj = c(1.04, -0.6))
))

edges <- list(
  list(text = "monocots", tips = "Oryza sativa"),
  list(text = "malvids", tips = "Arabidopsis thaliana")
)

whiteboard_edges <- c(edges, list(
  list(text = "fabids", tips = c("Tetraena mongolica", "Zygophyllum stapffii"))
))

dated_edges <- c(edges, list(
  list(text = "Zygophyllaceae",
       tips = c("Tetraena mongolica", "Zygophyllum stapffii"))
))

plot_dated_tree(whiteboard, whiteboard_nodes, whiteboard_edges,
                file.path(out_dir, "species_tree_whiteboard.png"))
plot_dated_tree(dated, dated_nodes, dated_edges,
                file.path(out_dir, "species_tree_dated.png"))
