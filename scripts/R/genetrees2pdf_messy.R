library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(readxl)
library(treeio)
library(magrittr)
library(RColorBrewer)


isolates = read_xlsx("/home/amsesk/data/pursuit/more_hmm/fasttree/../../Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, TAX.GROUP, LTP)

##### Batch print annotated trees #####
raw_path = "~/DATA/phylogeny/round1_raw/fast_gene_trees_renamed/"
filt_path = "/home/amsesk/data/pursuit/more_hmm/fasttree"

path = filt_path

files = list.files(path = path)

scores = read_delim("/home/amsesk/data/pursuit/more_hmm/unaln/hit_report_all.tsv", delim="\t", col_names=c("protein", "marker", "evalue", "score")) %>%
  separate("protein", into=c("LTP","protein"), sep="[|]") %>%
  left_join(isolates) %>%
  unite("tiplab", c("LTP","protein"), sep="|") %>%
  unite("tiplab", c("SPECIES.TREE.LABEL","tiplab"), sep="&") %>%
  mutate(tiplab_display = tiplab)

#scores$TAX.GROUP = plyr::revalue(scores$TAX.GROUP, c(
#    "chytrid" = "bold",
#    "non_chytrid_fungus" = "italic",
#    "outgroup" = "plain"
#  ))

#scores %<>%
#  mutate(TAX.GROUP = as.character(TAX.GROUP))

setwd("/home/amsesk/data/pursuit/more_hmm/fasttree/pdf")
for (f in files) {
  print(f)
  tree = read.newick( paste(path, f, sep="/") )

  tree = midpoint.root(tree)
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]tre", "", f)
  d = scores %>%
    filter(marker == this_marker)
  annotated_tree = ggtree(tree, branch.length = 0.1) %<+% d +
    ggtitle(this_marker) +
    #geom_tiplab(aes(label = paste(tiplab_display, score, sep = " - "), color=score, linetype = TAX.GROUP), cex=3) +
    geom_tiplab(aes(subset = TAX.GROUP == "chytrid", label = paste(tiplab_display, score, sep = " - "), color=score), size=2.5, fontface="bold")  +
    geom_tiplab(aes(subset = TAX.GROUP == "zoosporic", label = paste(tiplab_display, score, sep = " - "), color=score), size=2.5, fontface="bold.italic")  +
    geom_tiplab(aes(subset = TAX.GROUP == "non_chytrid_fungus", label = paste(tiplab_display, score, sep = " - "), color=score), size=2.5, fontface="italic")  +
    geom_tiplab(aes(subset = TAX.GROUP == "outgroup", label = paste(tiplab_display, score, sep = " - "), color=score), size=2.5)  +
    scale_color_continuous(type="viridis") +
    theme_tree()

  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, 10)

  ggsave(filename = paste(f, ".score_filter.scored", ".pdf", sep=""),
         plot = annotated_tree,
         device = "pdf",
         width=8.5,
         height=11,
         #width=17,
         #height=(14*((ntips)/140)),
         units="in",
         limitsize = FALSE)
}
annotated_tree

##### With node numbers #####

scores = read_delim("~/DATA/phylogeny/hit_report_all.csv", delim="\t", col_names=c("protein", "marker", "evalue", "score")) %>%
  separate("protein", into=c("LTP","protein"), sep="[|]") %>%
  left_join(isolates) %>%
  unite("tiplab", c("LTP","protein"), sep="|") %>%
  unite("tiplab", c("SPECIES.TREE.LABEL","tiplab"), sep="&") %>%
  mutate(tiplab_display = tiplab) %>%
  #mutate(tiplab_display = function(x) paste(strsplit(x, split="_")[[1]][1], strsplit(x, split="_")[[1]][2]))
  mutate(tiplab_display = str_extract(tiplab_display, "^[A-Za-z]+[_][A-Za-z]+"))

setwd("/home/aimzez/DATA/phylogeny/round5_final_score_filt/pdf_nodenum")
p=1
for (f in files) {
  print(f)
  tree = read.newick( paste(path, f, sep="/") )

  tree = midpoint.root(tree)
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)

  edge = data.frame(tree$edge, edge_num=1:nrow(tree$edge))
  colnames(edge)=c("parent", "node", "edge_num")

  d = scores %>%
    filter(marker == this_marker)

  annotated_tree = ggtree(tree, branch.length = 0.1) %<+% d %<+% edge +
    geom_label(aes(x=branch, label=edge_num), cex=2, label.padding=unit(0.05, "cm")) +
    geom_tiplab(aes(subset = TAX.GROUP == "chytrid", label = paste(tiplab_display, score, sep = " - "), color=score), cex=3, fontface="bold")  +
    geom_tiplab(aes(subset = TAX.GROUP == "zoosporic", label = paste(tiplab_display, score, sep = " - "), color=score), cex=3, fontface="bold.italic") +
    geom_tiplab(aes(subset = TAX.GROUP == "non_chytrid_fungus", label = paste(tiplab_display, score, sep = " - "), color=score), cex=3, fontface="italic")  +
    geom_tiplab(aes(subset = TAX.GROUP == "outgroup", label = paste(tiplab_display, score, sep = " - "), color=score), cex=3)  +
    scale_color_continuous(type="viridis") +
    ggtitle(label=this_marker)
    theme_tree()

  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, x_max*1.5)

  p=p+1
  if (p == 6) {
    break
  }
  #ggsave(filename = paste(f,".score_filter.scored.nodenumbers", ".pdf", sep=""), plot = annotated_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
}
annotated_tree

#############################################
##### Problematic taxon-annoated trees ######
#############################################

##### Without node numbers ####
setwd("/home/aimzez/DATA/phylogeny/round8_mancur_of_round7/pdf")
path = "/home/aimzez/DATA/phylogeny/round8_mancur_of_round7/fast_gene_trees_renamed"
files = list.files(path = path)
p=1
problematic_taxa = c("Mdap" = 1, "Falb" = 2, "Mbre" = 3, "Olpbor1" = 4)
tipcol = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP) %>%
  mutate(tip_color = ifelse(LTP %in% names(problematic_taxa), problematic_taxa[LTP], 0)) %>%
  mutate(tip_color = as.factor(tip_color))


myColors <- brewer.pal(5,"Set1")
myColors[1] = "#000000"
names(myColors) = levels(tipcol$tip_color)
scale = scale_colour_manual(name = "tip_color",values = myColors)

for (f in files) {
  print(f)
  tree = read.newick( paste(path, f, sep="/") )

  tree = midpoint.root(tree)
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)

  these_tipcol = as_tibble(tree$tip.label) %>%
    separate(value, into=c("SPECIES.TREE.LABEL", "pid"), sep="&") %>%
    left_join(tipcol) %>%
    unite(label, SPECIES.TREE.LABEL:pid, sep="&")

  annotated_tree = ggtree(tree, branch.length = 0.1) %<+% these_tipcol +
    geom_tiplab(aes(subset = tip_color == 0, color=tip_color), cex=3)  +
    geom_tiplab(aes(subset = tip_color != 0, color=tip_color), cex=3, fontface="bold")  +
    ggtitle(label=this_marker) +
    guides(
      color = FALSE
    ) +
    scale +
  theme_tree()

  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, x_max*1.5)
  #break
  ggsave(filename = paste(f,".round8.probtaxa_HL", ".pdf", sep=""), plot = annotated_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
 }
annotated_tree

##### WITH node numbers ####
setwd("/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/pdf_nodenum")
path = "/home/aimzez/DATA/phylogeny/round7_mancur_of_round5_spikein_rm_poly/fast_gene_trees_renamed"
files = list.files(path = path)
p=1
problematic_taxa = c("Mdap" = 1, "Falb" = 2, "Mbre" = 3, "Olpbor1" = 4)
tipcol = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP) %>%
  mutate(tip_color = ifelse(LTP %in% names(problematic_taxa), problematic_taxa[LTP], 0)) %>%
  mutate(tip_color = as.factor(tip_color))


myColors <- brewer.pal(5,"Set1")
myColors[1] = "#000000"
names(myColors) = levels(tipcol$tip_color)
scale = scale_colour_manual(name = "tip_color",values = myColors)

for (f in files) {
  print(f)
  tree = read.newick( paste(path, f, sep="/") )

  tree = midpoint.root(tree)
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)

  these_tipcol = tibble::as_tibble(tree$tip.label) %>%
    separate(value, into=c("SPECIES.TREE.LABEL", "pid"), sep="&") %>%
    left_join(tipcol) %>%
    unite(label, SPECIES.TREE.LABEL:pid, sep="&")

  annotated_tree = ggtree(tree, branch.length = 0.1) %<+% these_tipcol +
    geom_tiplab(aes(subset = tip_color == 0, color=tip_color), cex=3)  +
    geom_tiplab(aes(subset = tip_color != 0, color=tip_color), cex=3, fontface="bold")  +
    geom_text2(aes(subset=!isTip, label=node, fill=muted("red")), hjust=-.2, cex=2.6) +
    ggtitle(label=this_marker) +
    guides(
      color = FALSE
    ) +
    scale +
    theme_tree()

  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, x_max*1.5)
  #break
  ggsave(filename = paste(f,".round7.probtaxa_HL.nodenum", ".pdf", sep=""), plot = annotated_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
}
annotated_tree
