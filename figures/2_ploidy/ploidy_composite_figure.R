library(tidyverse)
library(scales)
library(tublerone)
library(phytools)
library(ggtree)
#library(ggimage)
library(readxl)
library(plotrix)
library(patchwork)

PATH_PREFIX = "/scratch/amsesk/pursuit/ploidy/"
SHEET_NAME = "../Pursuit_Isolates.xlsx"

weighted.sd = function(values, weights) {
  n_values = length(values)
  wm = weighted.mean(values, weights)
  val_min_wm_sq = (values - wm)^2
  the_sum = sum(val_min_wm_sq)
  v = the_sum / n_values
  sd = sqrt(v)
  return(sd)
}


colpal = c(
  "lightgrey",
  "red",
  "blue"
)

#### Kmer pane ####
setwd(file.path(PATH_PREFIX, "kmerhist"))
rambr = read_delim("Ramicandelaber_brevisporus_CBS_109374.Rambr1.v1_23.khist", delim="\t") %>%
  rename(Depth = `#Depth`)
spipal = read_delim("Spizellomyces_sp._palustris_CBS455.65_23.khist", delim="\t") %>%
  rename(Depth = `#Depth`)
rozal = read_delim("Rallo_23.khist", delim="\t") %>%
  rename(Depth = `#Depth`)

kmer_pane = ggplot(rambr) +
  geom_line(aes(x = Depth, y = logScale), color="blue", alpha=1.0) +
  geom_area(aes(x = Depth, y = logScale), fill="blue", alpha=0.4) +
  geom_line(data=spipal, aes(x=Depth, y=logScale), color="red", alpha=1.0) +
  geom_area(data=spipal, aes(x = Depth, y = logScale), fill="red", alpha=0.4) +
  geom_line(data=rozal, aes(x=Depth, y=logScale), color="darkgreen", alpha=1.0) +
  geom_area(data=rozal, aes(x = Depth, y = logScale), fill="darkgreen", alpha=0.4) +
  scale_x_continuous(expand = c(0,0), limits = c(25,400)) +
  scale_y_continuous(expand = c(0,0), limits = c(0,7.5e6)) +
  theme_bw() +
  ylab("kmer Counts") +
  xlab("kmer Depth") +
  theme (
    axis.title = element_text(size=8),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid = element_blank()
  )

#### AF pane ####
setwd(file.path(PATH_PREFIX, "snp_stats"))
cali_af = read_delim("California_12.snp_stats.tsv", delim="\t")
spipal_af = read_delim("Spizellomyces_sp._palustris_CBS455.65.snp_stats.tsv", delim="\t")
rozal_af = read_delim("Rallo.snp_stats.tsv", delim="\t")

af_pane = ggplot(cali_af) +
  geom_density(aes(x=p), fill=muted("blue"), alpha=0.4, color="blue") +
  geom_density(data=spipal_af, aes(x=p), fill=muted("red"), alpha=0.4, color= "red") +
  geom_density(data=rozal_af, aes(x=p), fill="darkgreen", alpha=0.4, color = "darkgreen") +
  scale_x_continuous(expand = c(0,0), limits = c(0,1)) +
  scale_y_continuous(expand = c(0,0)) +
  xlab("Reference Allele Frequency") +
  ylab("Density") +
  theme_bw() +
  theme (
    axis.title = element_text(size=8),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid = element_blank()
  )

#### Ploidy Scatter ####
l50_genome_sizes = read_delim(file.path(PATH_PREFIX, "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, l50_assembly_length)
isolates = read_xlsx(file.path(PATH_PREFIX, SHEET_NAME))
ploidy_df = read_delim(file.path(PATH_PREFIX, "Pursuit_Phylo_Traits.tsv"), delim="\t") %>%
  select(SPECIES.TREE.LABEL, coding) %>%
  rename(ploidy = coding)

snp_counts = read_delim(file.path(PATH_PREFIX, "all.snp_contig_counts_strainified.tsv"), delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  filter(ploidy_file_prefix == "Anasp1")

snp_densities = read_delim(file.path(PATH_PREFIX, "all.snp_contig_counts_strainified.tsv"), delim="\t", col_names = F) %>%
  rename(ploidy_file_prefix = X1, contig = X2, num_snps = X3, snp_density = X4, contig_length = X5) %>%
  group_by(ploidy_file_prefix) %>%
  summarise(snp_density_mean = weighted.mean(snp_density, contig_length), stdev = weighted.sd(snp_density, contig_length)) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  select(SPECIES.TREE.LABEL, snp_density_mean, stdev)

binom_expect = read_delim(file.path(PATH_PREFIX, "all_AF_ranges_from_binom.tsv"), delim="\t", col_names=F) %>%
  rename(ploidy_file_prefix = X1, mean_coverage = X2, p_in_binom_expect = X3) %>%
  left_join(isolates %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix)) %>%
  select(-ploidy_file_prefix) %>%
  select(SPECIES.TREE.LABEL, mean_coverage, p_in_binom_expect)

combined = binom_expect %>%
  left_join(snp_densities) %>%
  left_join(ploidy_df) %>%
  filter(mean_coverage > 10)

ploidy_scatter = ggplot(combined, aes(x = p_in_binom_expect, y = snp_density_mean, fill = ploidy)) +
  geom_hline(yintercept = 0, color="black", linetype="dashed") +
  geom_vline(xintercept=0, color="black", , linetype="dashed") +
  stat_ellipse(data=subset(combined, ploidy %in% c("2", "1")), geom="polygon", aes(fill=ploidy, color=ploidy), alpha=0.2, type="norm") +
  geom_errorbar(aes(ymax = snp_density_mean + stdev, ymin = snp_density_mean - stdev), alpha=0.30) +
  geom_point(pch = 21, cex = 2.5) +
  scale_fill_manual(values = colpal) +
  scale_color_manual(values = colpal[2:3]) +
  scale_x_continuous(expand = c(0,0.05)) +
  scale_y_continuous(expand = c(0,0)) +
  guides(fill = "none", color= "none") +
  xlab("SNPs in Expected") +
  ylab("Mean SNP Densitiy") +
  theme_bw() +
  theme(
    axis.title = element_text(size=8),
    plot.margin = unit(c(0,0,0,0), "in"),
    panel.grid = element_blank()
  )
#### Tree pane ####
isolates = read_delim(file.path(PATH_PREFIX, "Pursuit_Phylo_Traits.tsv"), delim = "\t") %>%
  select(SPECIES.TREE.LABEL, coding) %>%
  rename(ploidy = coding)

ploidy_tree = read.newick(file.path(PATH_PREFIX, "../combined_tree_filtered_support_RENAMED.tre"))

plt_ploidy_tree = ggtree(ploidy_tree, layout = "circular") %<+% isolates

ntips = length(ploidy_tree$tip.label)
nnodes = ploidy_tree$Nnode

#### Color EDGE by ploidy of descendant tips ####
#edgecol = matrix(nrow=ntips+nnodes, ncol=2)
#i=1
#for (n in seq(1, ntips+nnodes, 1)) {
#  edgecol[i,1] = n
#  desc = phytools::getDescendants(ploidy_tree, n)
#  desc_tips = desc[which(desc <= ntips)]
#  ploidy_of_desc_tips = sapply(desc_tips, function(x) plt_ploidy_tree$data[x,]$ploidy)
#  if (all(ploidy_of_desc_tips %in% c("?"))) {
#    edgecol[i,2] = "black"
#  } else {
#    if (all(ploidy_of_desc_tips %in% c("2.0", "?"))) {
#      edgecol[i,2] = "blue"
#    } else if (all(ploidy_of_desc_tips %in% c("1.0", "?"))) {
#      edgecol[i,2] = "red"
#    } else {
#      edgecol[i,2] = "purple"
#    }
#  }
#  i= i+1
#}

#edgecol = as_tibble(edgecol) %>%
#  rename(node = V1, edge_ploidy = V2) %>%
#  mutate(node = as.double(node))
#################################################

plt_ploidy_tree = plt_ploidy_tree +
  scale_fill_manual(values=colpal) +

  geom_hilight(node = 155, fill = "grey", alpha=0.6, extend = 100) +
  #geom_cladelabel(node = 155, label = "Chytridiomycota", fontsize = 6, offset=25, offset.text = 25) +

  #geom_hilight(node = 149, fill = "grey", alpha=0.4) +
  #geom_cladelabel(node = 149, label = "Neocallimastigomycota", fontsize = 6, offset=25, offset.text = 800) +

  geom_hilight(node = 153, fill = "grey", alpha=0.4, extend = 100) +
  #geom_cladelabel(node = 154, label = "monoblepharomycota", fontsize = 6, offset=25, offset.text = 800) +

  #geom_hilight(node = 241, fill = "grey", alpha=0.4, extend = 100) +
  #geom_cladelabel(node = 241, label = "Mucoromycota", fontsize = 6, offset=25, offset.text = 300) +

  geom_hilight(node = 229, fill = "grey", alpha=0.4, extend = 100) +
  #geom_cladelabel(node = 229, label = "Dikarya", fontsize = 6, offset=25, offset.text = 25) +

  #geom_hilight(node = 256, fill = "grey", alpha=0.4, extend = 100) +
  #geom_cladelabel(node = 256, label = "Zoopagomycota", fontsize = 6, offset=25, offset.text = 25) +

  geom_hilight(node = 266, fill = "grey", alpha=0.4, extend = 100) +
  #geom_cladelabel(node = 266, label = "Blastocladiomycota", fontsize = 6, offset=25, offset.text = 25) +

  #geom_hilight(node = 141, fill = "grey", alpha=0.4, extend=100) +
  #geom_cladelabel(node = 141, label = "Cryptomycota", fontsize = 6, offset=25, offset.text = 25)

  geom_hilight(node = 272, fill = "grey", alpha=0.4, extend=100) +
  #geom_cladelabel(node = 141, label = "Outgroup", fontsize = 6, offset=25, offset.text = 25)


  geom_tippoint(aes(fill=ploidy), color="black", pch=21, size=2)

plt_ploidy_tree + geom_text(aes(x=x, y=y, label=node))



#### Marginal Ancestral State Reconstruction ####
states = plt_ploidy_tree$data %>%
  select(label, isTip, ploidy) %>%
  filter(isTip) %>%
  select(-isTip) %>%
  mutate(ploidy = as.character(ploidy)) %>%
  mutate(ploidy = factor(ploidy, levels=c("1", "2")))
states_lst = states$ploidy
names(states_lst) = states$label

fitER = ace(states_lst, ploidy_tree, method = "ML", model="ARD", type="discrete")
anc_states_probs = as_tibble(fitER$lik.anc)

na_frame = cbind(rep(NA,137),rep(NA,137))
colnames(na_frame) = c("1", "2")

anc_states_final = rbind(na_frame,anc_states_probs) %>%
  mutate(node = seq(1, ntips+nnodes, 1)) %>%
  select(node, `1`, `2`)

plt_ploidy_tree$data = plt_ploidy_tree$data %>%
  left_join(anc_states_final)

marginal_tree_pane = plt_ploidy_tree +
  #aes(color=`2.0`) +
  geom_point(data = subset(plt_ploidy_tree$data, !isTip), aes(x=x, y=y, size=`2`), color = "blue", pch=16, alpha=0.3) +
  geom_point(data = subset(plt_ploidy_tree$data, !isTip), aes(x=x, y=y, size=`1`), color = "red", pch=16, alpha=0.3) +
  #scale_color_gradient2(low="red", mid = "purple", high="blue", midpoint=0.5) +
  scale_radius(range=c(0,4)) +
  guides(size=F, fill=F) +
  theme (
    plot.margin = unit(c(0,0,0,0), "cm"),
    plot.background = element_blank()
  )
marginal_tree_pane

### Pies on tree for supplemental
p = ggtree(ploidy_tree) %<+% isolates +
  geom_tiplab(size = 2.5, align = T, offset = 100, linesize=0.1, linetype="dashed") +
  xlim(0,1750) +
  geom_tippoint(aes(x=x+50, fill = ploidy), pch=21, size =2 ) +
  scale_fill_manual(values=colpal)
pies = nodepie(anc_states_final, cols=2:3, color = c("red", "blue"))
pies_on_tree = inset(p, pies, width =0.08, height =0.08)

pies_tree_cleaned_tiplabs = pies_on_tree
pies_tree_cleaned_tiplabs$data = pies_tree_cleaned_tiplabs$data %>%
  mutate(label = ifelse(isTip, gsub("[_.]v[0-9][.]*[0-9]*", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[.]LCG", "", label), label)) %>%
  mutate(label = ifelse(isTip, gsub("[_]", " ", label), label)) %>%
  mutate(label = gsub("LCG$", "", label)) %>%
  separate(label, c("genus", "species", "strain"), sep = " ", extra = "merge") %>%
  unite(col = "genus_species", genus, species, sep=" ", remove = F)

genus_species_duplicates = pies_tree_cleaned_tiplabs$data %>%
  filter(isTip) %>%
  select(genus_species) %>%
  group_by(genus_species) %>%
  summarise(occurances = n()) %>%
  filter(occurances > 1) %>%
  pull(genus_species)

pies_tree_cleaned_tiplabs$data = pies_tree_cleaned_tiplabs$data %>%
  mutate(label = ifelse(species == "sp.", paste(genus, species, strain), paste(genus, species))) %>%
  mutate(label = ifelse(genus_species %in% genus_species_duplicates, paste(genus, species, strain), label))


ggsave("/scratch/amsesk/Dropbox/Suppl_Ploidy_Tree_Pies.pdf",
       plot = pies_tree_cleaned_tiplabs,
       width = 8.5,
       height = 11,
       units = "in",
       device = cairo_pdf)

#### Marginal Ancestral State Reconstruction Important Node Bars
important = anc_states_final[c(140,155,149,241,229,266,258,153),]
important$group = c("Fungi", "Chytridiomycota", "Neocallimastigomycota", "Mucoromycota", "Dikarya", "Blastocladiomycota", "Zoopagomycota", "Monoblepharomycota")

melt = important %>%
  select(-node) %>%
  gather("key", "value", -group) %>%
  mutate(group = factor(group, levels = important$group[c(1,6,2,8,3,7,4,5)]))

marginal_anc_bars = melt %>%
  mutate(type = "marginal")

anc_state_bars = ggplot(marginal_anc_bars, aes(x= group, y=value, fill = key)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label = paste(round(value*100, 2), "%", sep="")), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values=c("red", "blue")) +
  theme_bw() +
  guides(fill=F) +
  ylab ("Percentage") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size=6)
  )

marginal_tree_pane + anc_state_bars

##################################################
##################################################
##################################################

#### MCMC Ancestral State Simulation
states = plt_ploidy_tree$data %>%
  select(label, isTip, ploidy) %>%
  filter(isTip) %>%
  select(-isTip) %>%
  mutate(ploidy = as.character(ploidy)) %>%
  mutate(ploidy = factor(ploidy, levels=c("1", "2")))

# Separate missing states from coded states
missing_states = states %>%
  filter(is.na(ploidy))
missing_labels = missing_states$label

coded_states = states %>%
  drop_na
coded_states_lst = coded_states$ploidy
names(coded_states_lst) = coded_states$label

# Make a tree that only has coded states to calculate simmap
coded_only_tree = ploidy_tree
coded_only_tree = drop.tip(coded_only_tree, tip = missing_states$label)
length(coded_only_tree$tip.label)

# Make the simmap based on coded states only
coded_states_simmap = make.simmap(coded_only_tree, coded_states_lst, model = "ARD")
Q = coded_states_simmap$Q

cols = setNames(c("red","blue"), c("1", "2"))

# Check to initial simmap tree
plot(coded_states_simmap, colors=cols, fsize=0.8, ftype="i")

# Generate matrix for known coded tips
x = coded_states %>%
  mutate(`1` = ifelse(ploidy == "1", 1, 0)) %>%
  mutate(`2` = ifelse(ploidy == "2", 1, 0)) %>%
  select(-ploidy)
rnames = x %>%
  pull(label)
x = x %>%
  select(-label)
X = as.matrix(x)
rownames(X) = rnames

# Make matrix for missing tips with 50-50 state probabilities
missing_X = matrix(nrow=length(missing_labels), ncol=2)
rownames(missing_X) = missing_labels
colnames(missing_X) = c("1", "2")
missing_X[,1] = 0.5
missing_X[,2] = 0.5

# Combine the coded and missing matrices
combined_states = rbind(X, missing_X)

# Now perform the stochastic mapping on the whole tree, missing and coded
mtrees = make.simmap(ploidy_tree, combined_states, model = "ARD", nsim=1000)

# Summarise the mtrees and get node probabilities
mtrees_sum = summary(mtrees)
node_probs = as_tibble(mtrees_sum$ace) %>%
  mutate(node = rownames(mtrees_sum$ace)) %>%
  mutate(node = as.numeric(node))


# Copy original ploidy_tree to stochastic one
stochastic_tree_pane = plt_ploidy_tree

# Remove marginal ancestral state node values and ploidy
stochastic_tree_pane$data = stochastic_tree_pane %>%
  select(-`1`, -`2`)

# Add stochastic ancestral state node values and ploidy
stochastic_tree_pane$data = stochastic_tree_pane$data %>%
  left_join(node_probs)

stochastic_tree_pane$data

# Plot ploidy_tree with the stochastic node labels
stochastic_tree_pane = stochastic_tree_pane +
  #aes(color=`2`) +
  geom_point(data = subset(stochastic_tree_pane$data, !isTip), aes(x=x, y=y, size=`2`), color = "blue", pch=16, alpha=0.3) +
  geom_point(data = subset(stochastic_tree_pane$data, !isTip), aes(x=x, y=y, size=`1`), color = "red", pch=16, alpha=0.3) +
  scale_radius(range=c(0,6)) +
  guides(size=F, fill=F) +
  theme (
    plot.margin = unit(c(0,0,0,0), "cm"),
    plot.background = element_blank()
  )
stochastic_tree_pane

#### Stochastic Ancestral State Reconstruction Important Node Bars
stochastic_important = node_probs %>%
  filter(node %in% c(140,155,149,241,229,266,258,153)) %>%
  arrange(node)
stochastic_important$group = c("Fungi",
                               "Neocallimastigomycota",
                               "Monoblepharomycota",
                               "Chytridiomycota",
                               "Dikarya",
                               "Mucoromycota",
                               "Zoopagomycota",
                               "Blastocladiomycota")

melt = stochastic_important %>%
  select(-node) %>%
  gather("key", "value", -group) %>%
  mutate(group = factor(group, levels = stochastic_important$group[c(1,8,4,3,2,7,6,5)]))

mcmc_anc_bars = melt %>%
  mutate(type="mcmc")

stochastic_anc_state_bars = ggplot(mcmc_anc_bars, aes(x= group, y=value, fill = key)) +
  geom_bar(stat="identity", color="black") +
  geom_text(aes(label = paste(round(value*100, 2), "%", sep="")), position = position_stack(vjust = 0.5), size = 3) +
  scale_fill_manual(values=c("red", "blue")) +
  theme_bw() +
  guides(fill=F) +
  ylab ("Percentage") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size=6)
  )

patchwork_design = "
144
244
344
344
355
"
fig = kmer_pane + af_pane + ploidy_scatter + stochastic_tree_pane + stochastic_anc_state_bars + plot_layout(design=patchwork_design)


### Combine anc state bars from marginal and mcmc ASR
combined_anc_bars = rbind(marginal_anc_bars, mcmc_anc_bars) %>%
  spread(key, value) %>%
  mutate(group = factor(group, levels=c("Fungi", "Blastocladiomycota", "Chytridiomycota",
                                        "Monoblepharomycota", "Neocallimastigomycota",
                                        "Zoopagomycota", "Mucoromycota", "Dikarya"))) %>%
  gather(key = "key", "value", -type, -group)
combined_anc_bars$group = plyr::revalue(combined_anc_bars$group, c(
  "Fungi" = 1,
  "Blastocladiomycota" = 2,
  "Chytridiomycota" = 3,
  "Monoblepharomycota" = 4,
  "Neocallimastigomycota" = 5,
  "Zoopagomycota" = 6,
  "Mucoromycota" = 7,
  "Dikarya" = 8
))
combined_anc_bars = combined_anc_bars %>%
  mutate(group = as.numeric(group))

barwidth = 0.45
combined_anc_bars_plot = ggplot(combined_anc_bars, aes(x= group, y=value, group = type, fill = key)) +
  geom_bar(data = subset(combined_anc_bars, type=="marginal"),
           aes(x= group, y=value, group = type, fill = key),
           stat="identity", color="black",position="stack", width=barwidth) +
  geom_text(data = subset(combined_anc_bars, type=="marginal"),
            aes(label = paste(format(round(value*100, 2), nsmall = 2), "%", sep=""), x = group),
            position = position_stack(vjust = 0.5), size = 1.7) +
  geom_bar(data = subset(combined_anc_bars, type=="mcmc"),
           aes(x= group + barwidth + 0.035, y=value, group = type, fill = key),
           stat="identity", color="black", position="stack", width=barwidth, alpha=0.7) +
  geom_text(data = subset(combined_anc_bars, type=="mcmc"),
            aes(label = paste(format(round(value*100, 2), nsmall = 2), "%", sep=""), x = group + barwidth + 0.035),
            position = position_stack(vjust = 0.5), size = 1.7) +
  scale_fill_manual(values=c("red", "blue")) +
  theme_bw() +
  guides(fill=F) +
  ylab ("Percentage") +
  theme(
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    plot.margin = unit(c(0,0,0,0), "cm"),
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank()
  )

patchwork_design = "
144
244
344
344
355
"
fig = kmer_pane +
  af_pane +
  ploidy_scatter +
  marginal_tree_pane +
  combined_anc_bars_plot +
  plot_layout(design=patchwork_design)

fig
ggsave(filename = file.path("/home/amsesk/Dropbox", "ploidy_combined_draft_091321_toIllustrator.pdf"),
       plot = fig,
       width = 8.5,
       height = 5.5,
       device = "pdf")

#### Phytools density map - not sure what exactly is going on with this one
cols = setNames(c("red","green"), c("1", "2"))
obj<-densityMap(mtrees,colors=c('red', 'blue'),lwd=4,outline=TRUE,type="fan",fsize = 0.3)
setMap(obj, colors=cols)
obj


#### Ploidy Scatter Summary Stats ####
View(combined %>%
  group_by(ploidy) %>%
  summarise(p_in_binom_expect_min = min(p_in_binom_expect),
            p_in_binom_expect_max = max(p_in_binom_expect),
            p_in_binom_expect_mean = mean(p_in_binom_expect),
            p_in_binom_expect_sd = sd(p_in_binom_expect),
            snp_density_mean_min = min(snp_density_mean),
            snp_density_mean_max = max(snp_density_mean),
            snp_density_mean_mean = mean(snp_density_mean),
            snp_density_mean_sd = sd(snp_density_mean))
)
