# Load all required libraries into environment
library(phytools)
library(ggtree) #Need to install, google this on bioconductor
library(ggimage) # install.packages("ggimage")
library(tidyverse) # install.packages("tidyverse")

# Read a tree, in newick format
tree = read.newick(file = "~/dev/rabern/ASTRAL_20201027.tre")

# Read in your csv
the_characters = as_tibble(read.table(file = "~/dev/rabern/ASTRAL_20201027_tree_morph.csv", sep = ",", fill = TRUE, header = TRUE)) %>%
  rename(SPECIES.TREE.LABEL = Taxon) #rename `Taxon` to `SPECIES.TREE.LABEL`

# Get you all the possible sub-lists in an object
names(tree)

# To get the node numbers
ggtree(tree) +
  geom_tiplab() +
  geom_text(aes(label=node))

# Extract a clade into a new tree by node number
# Change 187 to whatever clade your interested in, get that number from plotting the tree with node numbers
subtree = extract.clade(phy = tree, node = 187)
subtree$tip.label

# Creating a new dataframe (ie tibble) in order of the tree tip labels...
# ... to which we add our character data
the_characters_better = as_tibble(tree$tip.label) %>%
  rename(SPECIES.TREE.LABEL = value) %>%
  left_join(the_characters) %>% # Join tiplabel dataframe with character matrix; make sure both share the matching column
  mutate(M = as.factor(M)) %>% # Change `M` to whatever character
  filter(SPECIES.TREE.LABEL %in% subtree$tip.label) # Remove rows where the SPECIES.TREE.LABEL is not in `subtree`

# Creating a new ggtree object from `subtree` and joining `the_characters_better` with our_new_tree$data
our_new_tree = ggtree(subtree) %<+% the_characters_better

# Add geometries to our ggtree object that we created above
our_new_tree = our_new_tree +
  geom_tiplab() + # Tip labels
  geom_point(data = subset(our_new_tree$data, isTip == TRUE), aes(x=x, y=y, fill=M), pch=21, size=3) + #Node points on the tips, change `M` to whatever character
  theme_tree2() # Canned theme that has some tree stuff about it

# phytools::ace wants the states as a named list, where the names are the tiplabels and the values are the states
# Let's convert our tibble columns into that format
ugly_list = the_characters_better$M # creating a vector from the values in the `M` column, change `M` if need be
names(ugly_list) = the_characters_better$SPECIES.TREE.LABEL # Setting the names of our vector to the tiplabels from `the_characters_better`

# Running the Ancestral State Reconstruction
fitER = ace(ugly_list, subtree, model="ER", type="discrete")

# This is the element of `fitER` where the likelihoods are stored
fitER$lik.anc

# Convert matrix of likelihoods to tibble
anc_stats = as_tibble(fitER$lik.anc)

# Generate NAs to fill in tip nodes, add more `rep(NA,ntips)` for each additional state (ie factor levels)
ntips = length(subtree$tip.label)
na_df = as_tibble(cbind(
  rep(NA,ntips), rep(NA,ntips), rep(NA,ntips), rep(NA,ntips)
  ))

# So rbind doesn't complain, make sure that both matrices have the same names
names(na_df) = names(anc_stats)

# Blindly concatenate na_df and anc_stats, as long as they have the same number of columns and the same `names`
final_anc_stats = rbind(na_df, anc_stats) %>%
  mutate(node = seq(1,(2*ntips)-1,1)) # Add a new column, called "node", that is a sequence of numbers from 1 to (2*ntips)-1

# Generate a list of pie charts based on the final_anc_stats tibble
# Change the cols argument based on how many states there are (i.e., dim(final_anc_stats)[2]-1)
pie_list = nodepie(final_anc_stats, cols=1:4)

# Add the pie charts to `our_new_tree`
inset(our_new_tree, pie_list, width=0.02, height=0.1)
