library(castor)
library(tidyverse)
library(phytools)
library(ggtree)

CHYTRID_PHYLO = "/home/amsesk/dev/Chytrid-Phylogenomics/"

isolates = read_delim(file.path(CHYTRID_PHYLO, "spreadsheets", "Pursuit_Phylo_Traits.tsv"), delim = "\t") %>%
  select(SPECIES.TREE.LABEL, coding) %>%
  rename(ploidy = coding) %>%
  mutate(ploidy = as.character(ploidy)) %>%
  mutate(ploidy = factor(ploidy, levels=c("1", "2")))

ploidy_tree = read.newick("~/dev/Chytrid-Phylogenomics/figures/1_tree/combined_tree_filtered_support_RENAMED.tre")
ploidy_tree = tublerone::phylo_root(ploidy_tree, outgroup = c("Drosophila_melanogaster.v6", "Capsaspora_owczarzaki_ATCC_30864.v2", "Monosiga_brevicolis_MX1.v1"))

plt_ploidy_tree = ggtree(ploidy_tree) %<+% isolates

ntips = length(ploidy_tree$tip.label)
nnodes = ploidy_tree$Nnode

states = plt_ploidy_tree$data %>%
  select(label, isTip, ploidy) %>%
  filter(isTip) %>%
  select(-isTip)

states_vec = c()
for (tip in ploidy_tree$tip.label) {
  this_state = states %>% filter(label == tip) %>% mutate(ploidy = as.numeric(ploidy)) %>% pull(ploidy)
  states_vec = c(states_vec, this_state)
}

max_par = castor::hsp_max_parsimony(ploidy_tree, states_vec, transition_costs = "all_equal", weight_by_scenarios = FALSE)

max_par_probs = max_par$likelihoods
max_par_cost = max_par$total_cost

max_par_probs = as_tibble(max_par_probs) %>%
  mutate(node = seq(1, ntips+nnodes, 1)) %>%
  rename("1" = V1, "2" = V2) %>%
  select(node, everything()) %>%
  filter(node > ntips)

plt_ploidy_tree = plt_ploidy_tree +
  geom_tiplab(size = 2.0, align = T, offset =150) +
  geom_point(aes(x=x, y=y, fill = ploidy), pch=23, color= "black") +
  scale_fill_manual(values=c('red', 'blue'), na.value='grey') +
  scale_x_continuous(limits = c(0,2000),
                     expand=c(0,0)) +
  scale_y_continuous(expand=c(0, 0.6))

pies = nodepie(data = max_par_probs, cols = 2:3, color=c('red', 'blue'))

pie_tree_out = inset(plt_ploidy_tree, pies, width=0.1, height=0.1)

ggsave(filename = "/Users/aimzez/dev/Chytrid-Phylogenomics/figures_pnas_revisions/parsimony_acr/ploidy_parsimony_acr_pie_tree_WBS-False.pdf",
       plot = pie_tree_out,
       device = "pdf",
       width = 8.5,
       height = 11,
       unit = "in")

### WeightByScenario = TRUE
plt_ploidy_tree = ggtree(ploidy_tree) %<+% isolates

max_par = castor::hsp_max_parsimony(ploidy_tree, states_vec, transition_costs = "all_equal")
node_states = max.col(max_par$likelihoods)
max_par_probs = node_states

max_par_probs = as_tibble(max_par_probs) %>%
  mutate(node = seq(1, ntips+nnodes, 1)) %>%
  rename("state" = value) %>%
  select(node, everything()) %>%
  filter(node > ntips)

plt_ploidy_tree$data = plt_ploidy_tree$data %>%
  left_join(max_par_probs) %>%
  mutate(state = ifelse(isTip, ploidy, state)) %>%
  mutate(state = as.factor(state))

plt_ploidy_tree = plt_ploidy_tree +
  geom_tiplab(size = 2.0, align = T, offset =150) +
  geom_point(data = subset(plt_ploidy_tree$data, !isTip), aes(x=x, y=y, fill = state), pch=21, color= "black") +
  geom_point(data = subset(plt_ploidy_tree$data, isTip), aes(x=x, y=y, fill = ploidy), pch=21, color= "black") +
  scale_fill_manual(values=c('red', 'blue'), na.value='grey') +
  scale_x_continuous(limits = c(0,2000),
                     expand=c(0,0)) +
  scale_y_continuous(expand=c(0, 0.6))

ggsave(filename = "/Users/aimzez/dev/Chytrid-Phylogenomics/figures_pnas_revisions/parsimony_acr/ploidy_parsimony_acr_pie_tree_WBS-True-Decide.pdf",
       plot = plt_ploidy_tree,
       device = "pdf",
       width = 8.5,
       height = 11,
       unit = "in")

### Color the edge according to the direction of ploidy transition
iter = seq(1,dim(ploidy_tree$edge)[1],1)
h2d = 0
d2h = 0
edgecol = c("meh","meh")
for (i in iter) {
  left_state = plt_ploidy_tree$data %>% filter(node == ploidy_tree$edge[i,1]) %>% pull(state)
  right_state = plt_ploidy_tree$data %>% filter(node == ploidy_tree$edge[i,2]) %>% pull(state)
  if ( (is.na(left_state) || is.na(right_state)) ) {
    print(paste(ploidy_tree$edge[i,]))
    next
  }
  if (left_state == 1 && right_state == 2) {
    new_row = c(as.numeric(ploidy_tree$edge[i,2]), "haploid->diploid")
    h2d = h2d + 1
  } else if (left_state == 2 && right_state == 1) {
    new_row = c(as.numeric(ploidy_tree$edge[i,2]), "diploid->haploid")
    d2h = d2h + 1
  } else {
    new_row = c(as.numeric(ploidy_tree$edge[i,2]), "unchanged")
  }
  edgecol = rbind(edgecol, new_row)
}
edgecol = as_tibble(edgecol) %>%
  filter(V1 != "meh") %>%
  rename(node = V1, state_change = V2) %>%
  mutate(node =as.numeric(node))

h2d
d2h

plt_ploidy_tree = ggtree(ploidy_tree) %<+% isolates
plt_ploidy_tree$data = plt_ploidy_tree$data %>%
  left_join(max_par_probs) %>%
  left_join(edgecol) %>%
  mutate(state_change = as.factor(state_change)) %>%
  mutate(state = as.factor(state))
plt_ploidy_tree = plt_ploidy_tree +
  geom_tiplab(aes(color = state_change), size = 2.0, align = T, offset =150) +
  geom_point(data = subset(plt_ploidy_tree$data, !isTip), aes(x=x, y=y, fill = state), pch=21, color= "black") +
  geom_point(data = subset(plt_ploidy_tree$data, isTip), aes(x=x, y=y, fill = ploidy), pch=21, color= "black") +
  scale_fill_manual(values=c('red', 'blue'), na.value='grey') +
  scale_x_continuous(limits = c(0,2000),
                     expand=c(0,0)) +
  scale_y_continuous(expand=c(0, 0.6)) +
  scale_color_manual(values = c("blueviolet", "hotpink", "black"), na.value = "grey")

plt_ploidy_tree = plt_ploidy_tree + aes(color = state_change)
plt_ploidy_tree

ggsave(filename = file.path(CHYTRID_PHYLO, "figures_pnas_revisions/parsimony_asr/ploidy_parsimony_acr_pie_tree_WBS-True-Decide_Colored.pdf"),
       plot = plt_ploidy_tree,
       device = "pdf",
       width = 8.5,
       height = 11,
       unit = "in")
