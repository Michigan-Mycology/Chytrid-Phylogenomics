library(ape)
library(phytools)
library(tidyverse)
library(ggtree)
library(readxl)
library(treeio)
library(magrittr)
library(RColorBrewer)
library(scales)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, TAX.GROUP, LTP)

##### Batch print annotated trees #####
path = "/home/aimzez/DATA/phylogeny_trimfix/round4_mancur_round3/fast_gene_trees_renamed"

files = list.files(path = path)

was_it_spiked = read_delim("~/DATA/phylogeny_trimfix/round4_mancur_round3/was_it_spiked.tsv", delim="\t") %>%
  rename(LTP = X1)

probs = c(
  "Monosiga_brevicolis_MX1.v1",
  "Olpidium_bornovanus_UCB_F19785.Olpbor1",
  "Fonticula_alba_ATCC_38817.v2",
  "Mitosporidium_daphniae_UGP3"
)
setwd("/home/aimzez/DATA/phylogeny_trimfix/round4_mancur_round3/pdf")
for (f in files) {
  tree = read.newick( paste(path, f, sep="/") )
  
  tree = midpoint.root(tree)
  ntips = length(tree$tip.label)
  this_marker = gsub("[.]aa[.]tre[.]renamed", "", f)
  
  wtms = was_it_spiked %>% 
    select(LTP, !!sym(this_marker)) %>%
    rename(tips = !!sym(this_marker)) %>%
    mutate(tips = str_split(tips, ";")) %>%
    unnest(cols=c(tips)) %>%
    left_join(isolates %>% select(LTP, SPECIES.TREE.LABEL)) %>%
    drop_na %>%
    unite("combined", c(SPECIES.TREE.LABEL, tips), sep="&") %>%
    .$combined
  
  d = as_tibble(tree$tip.label) %>%
    mutate(tiplab = value) %>%
    separate(tiplab, sep="&", into=c("tiplab", "prot")) %>%
    mutate(spiked = ifelse(tiplab %in% probs, "PNS", "NP" )) %>%
    mutate(spiked = ifelse(value %in% wtms, "PS", spiked)) %>%
    mutate(spiked = factor(spiked, levels = c("NP", "PNS", "PS")))

  annotated_tree = ggtree(tree, branch.length = 0.1) %<+% d +
    geom_tiplab(aes(color=spiked), cex=3) + 
    ggtitle(this_marker) +
    scale_color_manual(values=c("NP"="black", "PNS"="green4", "PS"="deeppink1")) +
    #geom_text2(aes(subset=!isTip, label=node, fill=muted("red")), hjust=-.2, cex=2.6) +
    theme_tree()
  
  x_max = layer_scales(annotated_tree)$x$range$range[2]
  annotated_tree = annotated_tree +
    ggplot2::xlim(0, x_max*1.5)
  #break
  ggsave(filename = paste(f, "spikeHL", ".pdf", sep=""), plot = annotated_tree, device = "pdf", width=17, height=(14*((ntips)/140)), units="in", limitsize = FALSE)
}
annotated_tree
