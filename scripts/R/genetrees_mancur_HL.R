library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(readxl)
library(treeio)
library(magrittr)
library(RColorBrewer)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, TAX.GROUP, LTP)

##### Batch print annotated trees #####
raw_path = "~/DATA/phylogeny_trimfix/round1_raw/fast_gene_trees_renamed/"
filt_path = "/home/aimzez/DATA/phylogeny_trimfix/round3_scorefilt_25perc_gaps/fast_gene_trees_renamed"

path = filt_path

files = list.files(path = path)

scores = read_delim("~/DATA/phylogeny/hit_report_all.tsv", delim="\t", col_names=c("protein", "marker", "evalue", "score")) %>%
  separate("protein", into=c("LTP","protein"), sep="[|]") %>%
  left_join(isolates) %>%
  unite("tiplab", c("LTP","protein"), sep="|") %>%
  unite("tiplab", c("SPECIES.TREE.LABEL","tiplab"), sep="&") %>%
  mutate(tiplab_display = tiplab)

dropped_tips = read_delim("~/DATA/phylogeny_trimfix/round1_manucur_dropped_tips.tsv", delim="\t", col_names = c("marker", "tips"))[-1,]

setwd("/home/aimzez/DATA/phylogeny_trimfix/round3_scorefilt_25perc_gaps/gap25_mancur_round1_pdf")
for (f in files) {
  #print(f)
  tree = read.newick( paste(path, f, sep="/") )
  
  tree = midpoint.root(tree)
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)
  dt = dropped_tips %>%
    filter(marker == this_marker) %>%
    .$tips
  dt = unlist(strsplit(dt, split = ", "))
  
  annot = as_tibble(tree$tip.label) %>%
    rename(tip = value) %>%
    mutate(group = tip %in% dt)

  annotated_tree = ggtree(tree, branch.length = 0.1) %<+% annot +
    geom_tiplab(aes(color=group), cex=3) + 
    ggtitle(this_marker) +
    theme_tree()
  
  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, x_max*1.5)
  
  ggsave(filename = paste(f, "mancur_HL", ".pdf", sep=""), plot = annotated_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
}
annotated_tree
