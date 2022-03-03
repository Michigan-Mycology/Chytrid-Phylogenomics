library(tidyverse)
library(ggtree)
library(phytools)
library(scales)
library(readxl)

tree = read.newick("~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_all_support.tre")
tree$node.label
n_tips = length(tree$tip.label)

### Eliminate support values on nodes whose tips are all the same species, cause Tim suggested this to reduce noise in the figure
for (i in seq(n_tips+1, n_tips+(n_tips-1), 1)) {
  descendant_tips = c()
  descendant_nodes = phytools::getDescendants(tree, i)
  for (n in descendant_nodes) {
    if (n <= n_tips) {
      descendant_tips = c(descendant_tips, tree$tip.label[n])
    }
  }
  tip_species = as_tibble(descendant_tips) %>%
    separate(value, into=c("genus", "species")) %>%
    unite(taxon, genus:species, sep = "_") %>%
    .$taxon

  if (length(unique(tip_species)) == 1) {
    tree$node.label[i-n_tips] = "DONT_SHOW"
  }
}

# write.tree(tree, file = "/home/aimzez/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_filtered_support.tre")

plt_tree = ggtree(tree)

tree_data = as_tibble(plt_tree$data$label) %>%
  rename(label = value) %>%
  mutate(bootstrap = label) %>%
  separate("bootstrap", into = c("bootstrap", "qpic", "gcf", "astral"), sep = "/") %>%
  mutate_at(vars(bootstrap, qpic, gcf, astral), as.numeric) %>%
  mutate(gcf = gcf/100)

nodenum = ggtree(tree) +
  geom_tiplab(cex=2.5) +
  geom_text(aes(label=node))+
  theme_tree2()
nodenum

### Plot the tree

final = ggtree(tree) %>%
  collapse(node=229) %>%
  collapse(node=241) %>%
  collapse(node=272)
final = final %<+% tree_data

final = final +
  geom_text2(aes(subset=(node == 229)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45) +
  geom_text2(aes(subset=(node == 229)), label = "Dikarya", cex=3.0, vjust=0.4, hjust = -0.5) +
  geom_text2(aes(subset=(node == 241)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45) +
  geom_text2(aes(subset=(node == 241)), label = "Mucoromycota", cex=3.0, vjust=0.4, hjust = -0.25) +
  geom_text2(aes(subset=(node == 272)), cex=8, label=intToUtf8(9668), hjust =.2,vjust=.45) +
  geom_text2(aes(subset=(node == 272)), label = "Non-Fungal Eukaryota", cex=3.0, hjust =-0.2,vjust=.45)

x_max = max(final$data$x[!is.na(final$data$x)])
final = final +
  theme_tree2() +
  geom_tiplab(cex=2.5)
  #geom_point(data = subset(final$data, !isTip & !is.na(qpic)), aes(x=x-x_max, y=y, fill=qpic), color="black", pch=21, cex=2.0) +
  #geom_point(data = subset(final$data, !isTip & !is.na(gcf)), aes(x=x-x_max-20, y=y, fill=gcf), color="black", pch=22, cex=2.0) +
  #geom_point(data = subset(final$data, !isTip & label != "Root/"), aes(x=x-x_max-40, y=y, fill=astral), color="black", pch=23, cex=2.0) +
  #scale_fill_continuous(type="viridis", na.value="white")

horiz_just = 35
vertical_just = 1
nodelab_pos = as_tibble(final$data) %>%
  select(node, label, branch.length, x, y, isTip) %>%
  filter(!isTip) %>%
  mutate(repos_node_label = ifelse(branch.length < 35, TRUE, FALSE)) %>%
  mutate(node_x = x - x_max) %>%
  mutate(node_y = y) %>%
  mutate(nodelab_x = ifelse(repos_node_label, x - x_max - horiz_just, x - x_max)) %>%
  mutate(nodelab_y = ifelse(repos_node_label, y + vertical_just, y)) %>%
  select(-x, -y, -isTip, -label, -branch.length)
final$data = final$data %>%
  left_join(nodelab_pos, by="node")

final = final +
  geom_point(data = subset(final$data, !isTip & !is.na(qpic)), aes(x=nodelab_x, y=nodelab_y, fill=qpic), color="black", pch=21, cex=2.0) +
  geom_point(data = subset(final$data, !isTip & !is.na(gcf)), aes(x=nodelab_x-20, y=nodelab_y, fill=gcf), color="black", pch=22, cex=2.0) +

  # Show astral support values
  #geom_point(data = subset(final$data, !isTip & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW"), aes(x=nodelab_x-40, y=nodelab_y, fill=astral), color="black", pch=23, cex=2.0) +

  # Shown only where clade was not represented in astral tree
  geom_point(data = subset(final$data, !isTip & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW" & is.na(astral)), aes(x=nodelab_x-40, y=nodelab_y, fill=astral), color="black", pch=23, cex=2.0) +
  scale_fill_continuous(type="viridis", na.value="red")

final = revts(final) +
  geom_segment(data = subset(final$data, !isTip & repos_node_label & label != "NA/NA/NA/Root" & label != "/" & label != "DONT_SHOW"), aes(x=node_x, y=node_y, xend=nodelab_x, yend=nodelab_y-0.35), alpha=1.0) +
  #geom_text(data = subset(final$data, !isTip), aes(x=nodelab_x, y=nodelab_y, label=label), size=1.8, color=muted("blue")) +
  scale_x_continuous(limits = c(-1100,800), labels = abs, breaks = c(-1000, -750, -500, -250, 0)) +
  xlab("Millions of Years Ago (Mya)") +
  guides(
    fill = FALSE
  )

#final = revts(final) +
#  scale_x_continuous(limits = c(-1100,800), labels = abs, breaks = c(-1000, -750, -500, -250, 0)) +
#  xlab("Millions of Years Ago (Mya)") +
#  guides(
#    fill = FALSE
#  )

final

##### With showing all astral #####
#ggsave(plot = final, file = "~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_nodepoint_collpased.png", device = "png", width = 8.5, height = 11, units = "in")

#With showing only where astral differed
ggsave(plot = final, file = "~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_nodepoint_collpased_ASTRAL_mismatch_RED.png", device = "png", width = 8.5, height = 11, units = "in")

#### Map ploidy on to tips ####
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP, STRAIN, ploidy_file_prefix, "ploidy (12-11-20)") %>%
  rename(inferred_ploidy = "ploidy (12-11-20)") %>%
  mutate(inferred_ploidy = factor(inferred_ploidy, levels = c("2N", "1N", "2N?", "1N?", "?", "Not_Yet_Inferred"))) %>%
  select(SPECIES.TREE.LABEL, inferred_ploidy)

isolates
colpal = c(
  "blue",
  "red",
  muted("blue"),
  muted("red"),
  "yellow",
  "grey"
)
final_pp = final %<+% isolates
final_pp = final_pp +
  geom_point(data = subset(final_pp$data, node <= n_tips), aes(x=x, y=y, color = inferred_ploidy), pch = 20, cex =3) +
  scale_color_manual(values=colpal)

ggsave(plot = final_pp, file = "~/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_nodepoint_collpased_ASTRAL_mismatch_RED_ploidyDots.png", device = "png", width = 8.5, height = 11, units = "in")

######
