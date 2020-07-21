library(tidyverse)
library(readxl)

tims_markers = c("10168at4751",
                 "10531at4751",
                 "12258at4751",
                 "12812at4751",
                 "13224at4751",
                 "13605at4751",
                 "14295at4751",
                 "7545at4751",
                 "8406at4751",
                 "8818at4751")

cutoffs = read_delim("~/dev/domtbl2unaln/lib/odb10_scores_cutoff", delim="\t", col_names=FALSE) %>%
  rename(marker = X1, cutoff = X2) %>%
  filter(marker %in% tims_markers)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP)

domtbl_scores = read_delim("~/DATA/phylogeny/no_busco/hit_report_all.csv", delim="\t", col_names=c("protein", "marker", "evalue", "score")) %>%
  separate("protein", into=c("LTP","protein"), sep="[|]") %>%
  left_join(isolates) %>%
  unite("tiplab", c("LTP","protein"), sep="|") %>%
  unite("tiplab", c("SPECIES.TREE.LABEL","tiplab"), sep="&") %>%
  filter(marker %in% tims_markers)

myplots = list()
for (t in tims_markers) {
  treef = paste(t,".aa.tre.renamed",sep="")
  tree = read.newick(paste("~/DATA/phylogeny/round3_scorefilt_recalc_monogenus/fast_gene_trees_filtered_renamed/", treef, sep=""))
  these_marker_scores = domtbl_scores %>%
    filter(marker == t)
  mitospor = these_marker_scores %>% 
    filter(str_detect(tiplab, "Mitosporidium")) %>%
    .$score
  print(mitospor)
  the_plot = ggplot(these_marker_scores, aes(x=score)) +
    geom_histogram(bins = 50, color="black", fill="lightgreen", alpha=0.8) +
    geom_vline(xintercept = cutoffs[cutoffs$marker == t,]$cutoff, color="red") + 
    geom_vline(xintercept = mitospor, color="blue") +
    ggtitle(t) +
    theme_bw()
  
  myplots[[t]] = the_plot
}

plt_pane(myplots, 5,2)
