library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(scales)

args = commandArgs(trailingOnly = TRUE)
tree_path = args[1]

#tree_path = "~/DATA/phylogeny_final/concat/treefiles"
#setwd("~/DATA/phylogeny_final/concat/treefiles/")

tree_files = list.files(path = tree_path)

for (f in tree_files) {
#f = "part_concat_UFB.treefile.renamed"
  print(f)
  tree = read.newick( paste(tree_path, f, sep="/") )
  
  ### Root tree, pick the node if not midpoint rooting
  tree = midpoint.root(tree)
  #tree = ape::root(tree, node = 226)
  #tree = ape::root(tree, node = 138)
  
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]tre[.]renamed", "", f)
  annotated_tree = ggtree(tree, branch.length = 0.1)
  
  ### Fix bootstrap value not being transcribed in reroot/midpoint.root
  annotated_tree$data$label = annotated_tree$data %>% mutate(label = ifelse(label == "", "100", label)) %>% .$label
  
  annotated_tree = annotated_tree +
    geom_tiplab(cex=3) +
    #geom_nodelab(cex=3) +
    #geom_nodelab(aes(label=node)) +
    #geom_nodepoint(aes(fill=as.numeric(label), cex=as.numeric(label)), alpha=0.5, shape=21, color="black") +
    #geom_nodepoint(aes(subset = label == 100, cex=as.numeric(label)), fill=muted("blue"), alpha=0.5, shape=21, color="black", show.legend = FALSE) +
    #geom_nodepoint(aes(subset = as.numeric(label) < 100 & as.numeric(label) > 0, cex=as.numeric(label)), fill=muted("red"), alpha=0.5, shape=21, color="black", show.legend = FALSE) +
    scale_fill_continuous(type="viridis", name="Bootstrap", limits=c(0,100)) +
    theme_tree() +
    labs(cex = "Bootstrap")
  
  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, x_max*1.5)

  ggsave(filename = paste(f, ".pdf", sep=""), plot = annotated_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
}
annotated_tree
