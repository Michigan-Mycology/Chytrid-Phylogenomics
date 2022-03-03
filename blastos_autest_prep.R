library(phytools)
library(ggtree)

best_tree = read.newick("/home/aimzez/DATA/pursuit/phylogeny_final/tims_trees/combined_tree_filtered_support.tre")
best_tree = ape::root(best_tree, outgroup=c("Monosiga_brevicolis_MX1.v1", "Drosophila_melanogaster.v6", "Capsaspora_owczarzaki_ATCC_30864.v2"))
ggtree(best_tree) +
  geom_nodelab(aes(label=node)) +
  geom_tiplab(cex=2.5) +
  xlim(0,1500)

terrestrial_fungi = 0

blasto_tips = best_tree$tip.label[getDescendants(best_tree, 266)[which(getDescendants(best_tree, 266) <= length(best_tree$tip.label))]]
best_tree_minus_blastos = drop.clade(best_tree, blasto_tips)
best_tree_minus_blastos = drop.tip(best_tree_minus_blastos, "NA")
blasto_clade = extract.clade(best_tree, 266)
blasto_clade$edge.length = rep(1, length(blasto_clade$edge.length))

ggtree(best_tree_minus_blastos) +
  geom_nodelab(aes(label=node)) +
  geom_tiplab(cex=2.5) +
  xlim(0,1500)

chy_aph =  bind.tree(best_tree_minus_blastos, blasto_clade, 139)
ggtree(chy_aph) +
  geom_nodelab(aes(label=node)) +
  geom_tiplab(cex=2.5)

#insert_loc = which(best_tree_minus_blastos$edge[,1] == 138 & best_tree_minus_blastos$edge[,2] == 139)
#chy_aph_edge = matrix(nrow = length(best_tree_minus_blastos$edge[,1]) + 1, ncol =2)

#chy_aph_edge[-insert_loc,] = best_tree_minus_blastos$edge
#chy_aph_edge[insert_loc+1,] = c(138,260)
#chy_aph_edge[insert_loc,] = c(260,139)
#chy_aph_edge

#chy_aph = best_tree_minus_blastos
#chy_aph$edge.length = rep(1, 258)
#chy_aph$edge = chy_aph_edge
#chy_aph$node.label[130] = "NA"
#chy_aph$Nnode = chy_aph$Nnode + 1
#chy_aph$edge.length[259] = 1

#chy_aph = bind.tip(chy_aph, "NewTip", edge.length = 1, where=260)
#chy_aph = bind.tree(chy_aph, blasto_clade, where = which(chy_aph$tip.label == "NewTip"))

chy_aph = tublerone::insert_clade(blasto_clade, best_tree_minus_blastos, node = 138, one_later_diverging_tip = "Rhizopus_delemar_RA_99-880")
after_olp = tublerone::insert_clade(blasto_clade, best_tree_minus_blastos, node = 218, one_later_diverging_tip = "Yarrow_lipolytica_CLIB122")
sis_chy = tublerone::insert_clade(blasto_clade, best_tree_minus_blastos, node = 139, one_later_diverging_tip = "Piromyces_sp._E2_v1.0")
ggtree(sis_chy) +
  geom_tiplab(cex=2.5) +
  geom_nodelab(aes(label=node))

write.tree(chy_aph, "~/DATA/pursuit/phylogeny_final/blasto_autest/chyaph.tre")
write.tree(after_olp, "~/DATA/pursuit/phylogeny_final/blasto_autest/afterolp.tre")
write.tree(sis_chy, "~/DATA/pursuit/phylogeny_final/blasto_autest/sischy.tre")
