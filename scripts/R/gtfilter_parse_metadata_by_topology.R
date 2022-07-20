library(tidyverse)

meta = read_delim("~/DATA/bonitomes/all.marker_orthodb_metadata.tsv", delim ="\t", col_names = F) %>%
  rename(marker = X1, annot_class = X2, accession = X3, desc = X4)

m2t = read_delim("~/DATA/bonitomes/gtfilter_marker2topology.tsv", delim = "\t", col_names = F) %>%
  rename(marker = X1, topology= X2)

comb = meta %>%
  left_join(m2t)
comb %>%
  filter(topology == "((0, 3), (1, 2))") %>%
  filter(annot_class == "OrthoDB") %>%
  group_by(desc) %>%
  summarize(n=n()) %>%
  arrange(desc(n))

comb %>%
  filter(topology == "(0, (3, (1, 2)))") %>%
  filter(annot_class == "OrthoDB") %>%
  group_by(desc) %>%
  summarize(n=n()) %>%
  arrange(desc(n))
