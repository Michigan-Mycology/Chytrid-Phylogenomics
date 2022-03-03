library(tidyverse)

scores = read_delim("~/DATA/phylogeny/hit_report_all.csv", delim="\t", col_names = FALSE) %>%
  mutate(X3 = as.double(X3))
colnames(scores) = c("protein", "marker", "evalue", "score")

alnlen = read_delim("~/DATA/phylogeny/aln_lengths.csv", delim="\t", col_names = FALSE)
colnames(alnlen) = c("marker", "protein", "alnlen")
alnlen %>%
  select(protein, marker, alnlen)

combined = scores %>%
  left_join(alnlen)
  #filter(marker %in% unique(scores$marker)[1:5])

ggplot(combined, aes(x=alnlen, y=score, color=marker, pch=marker)) +
  geom_point(cex=1.5, alpha=0.8)

ggplot(combined, aes(x=alnlen, y=score)) +
  geom_point(cex=0.5, alpha=0.8)

ggplot(combined, aes(x=score, y=log(evalue))) +
  geom_point(cex=0.5)
