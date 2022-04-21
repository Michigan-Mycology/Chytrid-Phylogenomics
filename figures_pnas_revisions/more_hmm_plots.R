library(tidyverse)
library(patchwork)
library(RColorBrewer)
library(readxl)
library(tublerone)
library(phytools)
library(ggtree)

d = read_delim("~/DATA/pnas_rev/more_hmms/hit_report_best.tsv", delim = "\t", col_names = F) %>%
  rename(prot = X1, gene = X2, evalue = X3, bitscore = X4)


### Bitscore Histograms
p=1
plts = list()
pal = brewer.pal(7, "Paired")
for (i in unique(d$gene)) {

  plts[[p]] = ggplot(data=subset(d, d$gene == i), aes(x = bitscore, fill = gene)) +
    geom_density(alpha=0.5) +
    scale_fill_manual(values = c(pal[p])) +
    theme_bw()

  p = p+1
}
bithists = plts[[1]] + plts[[2]] + plts[[3]] + plts[[4]] + plts[[5]] + plts[[6]] + plts[[7]] + plot_layout(ncol=2, nrow=4)

ggsave(filename = "~/DATA/pnas_rev/more_hmms/bitscore_hists.pdf",
       device = cairo_pdf,
       plot = bithists,
       width =8.5,
       height =11,
       unit = "in")

### Occupancy and bitscore heatmap

rename_heatmap_rows = function(heatmap, nameswitch) {
  new_names = c()
  i = 1
  for (n in heatmap$LTP) {
      new_names = c(new_names, nameswitch[n][[1]])
    i = i +1
  }
  return(new_names)
}

CHYTRID_PHYLO="/home/aimzez/work/Chytrid-Phylogenomics"

tree = read.newick(file.path(CHYTRID_PHYLO, "figures/1_tree", "combined_tree_filtered_support_RENAMED.tre"))
sheet = read_xlsx("~/DATA/pnas_rev/Pursuit_Isolates.xlsx")
old_to_new = sheet$SPECIES.TREE.LABEL
names(old_to_new) = sheet$LTP

hmd = d %>%
  separate(prot, c("LTP", "prot"), sep = "[|]") %>%
  group_by(gene) %>%
  mutate(bitscore_max = max(bitscore)) %>%
  ungroup() %>%
  mutate(bitscore_scaled = bitscore/bitscore_max)
hmd = hmd %>%
  mutate(label = rename_heatmap_rows(hmd, old_to_new))

plttree = ggtree(tree) +
  geom_tiplab(size=2) +
  xlim(0,1350)

hmplt = tublerone::ggtreeplot(ggtree = plttree, data = hmd, mapping = aes(x=gene, fill = bitscore_scaled)) +
  geom_tile() +
  scale_fill_continuous(type = 'viridis') +
  theme_bw() +
  coord_equal() +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.text.x = element_text(angle = 60, hjust=1, size =6),
    axis.ticks.y = element_blank(),
    axis.ticks.length.y = unit(0, "cm")
  )

plttree + hmplt

ggsave(filename = "~/DATA/pnas_rev/more_hmms/best_hits_unfiltered_tree_heatmap.pdf",
       device = cairo_pdf,
       plot = plttree+hmplt,
       width =8.5,
       height =11,
       unit = "in")

