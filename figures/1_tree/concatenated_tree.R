library(tidyverse)
library(phytools)
library(ggtree)
library(tublerone)
library(readxl)
library(RColorBrewer)
library(grid)
library(gridExtra)
library(scales)

treepath = "~/DATA/pursuit/phylogeny_final/concat/treefiles/unpart_NPB_cons_outgrouped.treefile"

tree = read.newick(treepath)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", sheet=1)

tree = root(tree, outgroup=c("Mbre","Cowc","Dmel","Falb"), resolve.root = TRUE)

old = as.list(isolates$LTP)
names(old) = isolates$LTP
old_to_new = lapply(old, FUN = function(x) isolates %>% select(SPECIES.TREE.LABEL, LTP) %>% filter (LTP == x) %>% .$SPECIES.TREE.LABEL)

orders = isolates %>%
  select(SPECIES.TREE.LABEL, ORDER)

tree = tublerone::phylo_rename_tips(tree, old_to_new)
write.tree(tree, file = "~/DATA/pursuit/phylogeny_final/concat/treefiles/unpart_concat_NFB.treefile.renamed")

phylum_nodes = read_delim("~/work/Chytrid-Phylogenomics/figures/1_tree/phylum_nodemap", delim="\t", col_names=FALSE) %>%
  rename(phylum=X1, node=X2, color=X3)
order_nodes = read_delim("~/work/Chytrid-Phylogenomics/figures/1_tree/order_nodemap", delim="\t", col_names=FALSE) %>%
  rename(order=X1, node=X2, color=X3)

ultra = read_xlsx("~/work/pursuit/sheets/Phylogenomics_Ultrastructure.xlsx") %>%
  mutate(Taxon = gsub("[ ]","_", Taxon))
ultra_merger = isolates %>%
  select(SPECIES.TREE.LABEL) %>%
  mutate(hold = SPECIES.TREE.LABEL) %>%
  separate(SPECIES.TREE.LABEL, into=c("genus", "species"), sep="[_]") %>%
  mutate(genus_merge = genus) %>%
  mutate(species_merge = species) %>%
  unite(gs_combined, c(genus_merge,species_merge), sep="_")
ultra_merger


rename_ultra = function(gs_combined, genus, ultra) {
  is_gs = gs_combined %in% ultra$Taxon
  if (is_gs) {
    ultraTax = ultra %>%
      filter(Taxon == gs_combined) %>%
      .$Taxon

    return(ultraTax)

  } else {
    is_g = genus %in% ultra$Taxon
    if(is_g) {
      ultraTax = ultra %>%
        filter(Taxon == genus) %>%
        .$Taxon

      return(ultraTax)
    } else {
      return("NO")
    }
  }
}

ultra_merged = ultra_merger %>%
  rowwise() %>%
  mutate(Taxon = rename_ultra(gs_combined, genus, ultra)) %>%
  select(hold, Taxon) %>%
  rename(SPECIES.TREE.LABEL = hold) %>%
  left_join(ultra, by="Taxon") %>%
  select(-Taxon)

### Figure out which characters are non-variable and exclude them
for (c in colnames(ultra_merged)[-1]) {
  vec = ultra_merged[,c] %>%
    na.omit
  if ( dim(unique(vec))[1] == 1) {
    ultra_merged = ultra_merged %>%
      select(-!!sym(c))
  }
}

ultra_merged_long = ultra_merged%>%
  gather(key="char", value="state", -SPECIES.TREE.LABEL)

ultra_merged_hm = as.data.frame(ultra_merged)
rownames(ultra_merged_hm) = ultra_merged_hm$SPECIES.TREE.LABEL
ultra_merged_hm = ultra_merged_hm[,-1]

ultra_merged

baseplt = ggplot()
i = 0
for (c in unique(ultra_merged$character)) {
  sub = subset(ultra_merged, character == c)
  pal = brewer.pal(n=dim(sub)[1], name="RdYlBu")
  baseplt = baseplt +
    geom_tile(data = sub, aes(x=character, y=SPECIES.TREE.LABEL, fill=state)) +
    scale_fill_manual(values=pal)
}
baseplt

#### Plot the Tree ####
plt_tree = ggtree(tree, size = 1.5) +
  #geom_text(aes(label=node), color="red", cex =3) +
  ggplot2::xlim(0, 4)

plt_tree$data = plt_tree$data %>%
  mutate(label = ifelse(label == "Root", NA, label)) %>%
  mutate(label = ifelse(label == "", NA, label))

for (i in seq(1, length(phylum_nodes$phylum), 1)) {
  plt_tree = plt_tree +
    geom_hilight(
      node = phylum_nodes[i,]$node,
      extendto = 4,
      fill = phylum_nodes[i,]$color,
      alpha=0.4)
    #geom_cladelabel(
      #node = phylum_nodes[i,]$node,
      #label = phylum_nodes[i,]$phylum,
      #color = phylum_nodes[i,]$color,
      #align=T,
      #offset=1.5)
}
plt_tree = plt_tree +
  geom_tiplab(cex=2.5, align=T) +
  geom_nodepoint(aes(size=as.numeric(label))) +
  scale_size_continuous(range=c(1.5,1.5)) +
  guides(size=F)

for (i in seq(1, length(order_nodes$order), 1)) {
  plt_tree = plt_tree +
    geom_cladelabel(
      node = order_nodes[i,]$node,
      label = order_nodes[i,]$order,
      color = order_nodes[i,]$color,
      align=T,
      offset=1.5)
}

plt_tree

ultra_merged_long = ultra_merged_long %>%
  mutate(char = as.factor(char)) %>%
  mutate(state = as.factor(state))

pts = plt_tree
for (c in colnames(ultra_merged)[-1]) {
  print(c)
  pts = pts +
    geom_facet(panel = c, data = subset(ultra_merged_long, char == c), geom = geom_point, mapping=aes(x = 1, color = state))
}
pts
#pts +
#  geom_facet(panel = "M", data = subset(ultra_merged_long, char == "M"), geom = geom_point, mapping=aes(x = 1, color = state)) +
#  geom_facet(panel = "N", data = subset(ultra_merged_long, char == "N"), geom = geom_point, mapping=aes(x = 1, color = state))

pts = pts +
  scale_color_manual(values = brewer.pal(n=11, name = "Set3"), na.value="white") +
  theme (
    panel.spacing.x = unit(0, "in")
  )

gp = ggplotGrob(pts)
facet.columns <- gp$layout$l[grepl("panel", gp$layout$name)]
gp$widths[facet.columns[1]] = unit(5, "in")
gp$widths[facet.columns[-1]] = rep(unit(0.09722, "in"), 36)
grid::grid.newpage()
grid::grid.draw(gp)


pal = brewer.pal(n=11, name="RdYlBu")
gheatmap(plt_tree, ultra_merged_hm, offset=0.5, legend_title = "Character State", color=NULL) +
  scale_fill_manual(values=pal)

ggplot(ultra_merged_long) +
  geom_point(data = subset(ultra_merged_long, char == "M"), aes(x=1, y=seq(1,140,1), color=state)) +
  geom_point(data = subset(ultra_merged_long, char == "N"), aes(x=2, y=seq(1,140,1), color=state))

#### Investigating some of Tim's trees and where the pertinent data is
library(treeio)
qs_tree_qpic_only = function (nodelabel) {
  spl = strsplit(nodelabel, split = "/")
  qpic = gsub(pattern = "qp[-]ic[=]", replacement = "", spl[[1]][2])
  qpic = format(round(as.double(qpic), 2), nsmall = 2)
  return (paste(spl[[1]][1], qpic, sep="/"))
}
concord_tree_gcf_only = function (nodelabel) {
  spl = strsplit(nodelabel, split = "/")
  return (spl[[1]][2])
}

rates_tree = read.newick("~/DATA/pursuit/phylogeny_final/tims_trees/best_r8s_tree_Run1.txt")
concord_tree = read.newick("~/DATA/pursuit/phylogeny_final/tims_trees/concord_chytrid_iqtree.cf.tree")
qs_tree = read.newick("~/DATA/pursuit/phylogeny_final/tims_trees/QS_iqtree_tree_renamed_reformat.tre")

#Rename tips of rates_tree so they match with other two trees
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", sheet=1)
old = as.list(isolates$LTP)
names(old) = isolates$LTP
old_to_new = lapply(old, FUN = function(x) isolates %>% select(SPECIES.TREE.LABEL, LTP) %>% filter (LTP == x) %>% .$SPECIES.TREE.LABEL)
rates_tree = tublerone::phylo_rename_tips(rates_tree, old_to_new)

# Only include bootstrap and qp-ic from qs_tree nodelabels
qs_tree$node.label = unlist(lapply(qs_tree$node.label, qs_tree_qpic_only))

#Add gCF (gene Concordance Factor), the second one, support values from concord_tree
concord_tree$node.label = unlist(lapply(concord_tree$node.label, concord_tree_gcf_only))

combined_tree = qs_tree
combined_tree$node.label = paste(combined_tree$node.label, concord_tree$node.label, sep="/")
combined_tree = root(combined_tree, resolve.root = TRUE, edgelabel = TRUE, outgroup = c(
  "Monosiga_brevicolis_MX1.v1",
  "Drosophila_melanogaster.v6",
  "Capsaspora_owczarzaki_ATCC_30864.v2"
))

# Make sure that the architectures of rates_tree and combined_tree are the same
comparePhylo(combined_tree, rates_tree, plot = T)

# All this to move the node labels from combined_tree to rates_tree (since the nodes/tips are numbered differently)
tips_and_support_by_node = function(tree) {
  n_tips = length(tree$tip.label)
  n_internal = tree$Nnode
  nodes_and_tips = n_tips + n_internal
  nodes_and_their_tips = list()
  for (internal_node_number in seq(n_tips+1, nodes_and_tips, 1)) {
    descendant_tips = c()
    descendant_nodes = phytools::getDescendants(tree, internal_node_number)
    for (n in descendant_nodes) {
      if (n <= n_tips) {
        descendant_tips = c(descendant_tips, tree$tip.label[n])
      }
    }
    nodes_and_their_tips[[as.character(internal_node_number)]] = list( tips = descendant_tips, support = tree$node.label[internal_node_number-n_tips] )
  }
  nodes_and_their_tips
}
the_map = tips_and_support_by_node(combined_tree)

astral_tree = read.newick("~/DATA/pursuit/astral_tree/ASTRAL_20201027.noBrlen.tre.renamed")
astral_tree = root(astral_tree, resolve.root = TRUE, edgelabel = TRUE, outgroup = c(
  "Monosiga_brevicolis_MX1.v1",
  "Drosophila_melanogaster.v6",
  "Capsaspora_owczarzaki_ATCC_30864.v2"
))
astral_map = tips_and_support_by_node(astral_tree)

for (i in seq(138, 273, 1)) {
  descendant_tips = c()
  descendant_nodes = phytools::getDescendants(rates_tree, i)
  for (n in descendant_nodes) {
    if (n <= 137) {
      descendant_tips = c(descendant_tips, rates_tree$tip.label[n])
    }
  }
  print (descendant_tips)

  for (p in names(the_map)) {
    if (length(descendant_tips) == length(the_map[[p]]$tips) & all(descendant_tips %in% the_map[[p]]$tips)) {
      rates_tree$node.label[i-137] = the_map[[p]]$support
    }
  }
  for (p in names(astral_map)) {
    if (length(descendant_tips) == length(astral_map[[p]]$tips) & all(descendant_tips %in% astral_map[[p]]$tips)) {
      rates_tree$node.label[i-137] = paste(rates_tree$node.label[i-137], astral_map[[p]]$support, sep="/")
    }
  }
  support_spl = unlist(strsplit(rates_tree$node.label[i-137], split="/"))
  if ( length(support_spl) == 3 ) {
    rates_tree$node.label[i-137] = paste(rates_tree$node.label[i-137], "NA", sep="/")
  }

  if ( length(support_spl) == 2 ) {
    rates_tree$node.label[i-137] = paste("NA/NA/NA", astral_map[[i-137]]$support, sep="/")
  }

}

### Write out `rates_tree`, which is really the combined tree with all support values, to newick
write.tree(phy = rates_tree, file="~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_all_support.tre")

# Plot the finalized, combined tree
horiz_just = 75
vertical_just = 1

final_tree = ggtree(rates_tree)  +
  geom_tiplab(cex=2.5) +
  #geom_nodelab(size=3, hjust=horiz_just, vjust=vertical_just) +
  xlab("Millions of Years Ago (Mya)") +
  theme_tree2()

nodelab_pos = as_tibble(final_tree$data) %>%
  select(node, label, branch.length, x, y, isTip) %>%
  filter(!isTip) %>%
  mutate(repos_node_label = ifelse(branch.length < 35, TRUE, FALSE)) %>%
  mutate(node_x = x - max(x)) %>%
  mutate(node_y = y) %>%
  mutate(nodelab_x = ifelse(repos_node_label, x - max(x) - horiz_just, x - max(x) - 50)) %>%
  mutate(nodelab_y = ifelse(repos_node_label, y + vertical_just, y + 0.75)) %>%
  select(-x, -y, -isTip, -label, -branch.length)
final_tree$data = final_tree$data %>%
  left_join(nodelab_pos, by="node")

final_tree = revts(final_tree) +
  geom_segment(data = subset(final_tree$data, !isTip & repos_node_label), aes(x=node_x, y=node_y, xend=nodelab_x, yend=nodelab_y-0.35), alpha=1.0, size=0.1) +
  geom_text(data = subset(final_tree$data, !isTip), aes(x=nodelab_x, y=nodelab_y, label=label), size=1.8, color=muted("blue")) +
  scale_x_continuous(limits = c(-1100,500), labels = abs)

ggsave(plot = final_tree, file = "~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree.pdf", device = "pdf", width = 8.5, height = 11, units = "in")
