library(tidyverse)
library(readxl)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, TAX.GROUP, LTP)
repr = read_xlsx("~/DATA/phylogeny/stepwise_taxon_repr.xlsx")

repr %<>%
  left_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything())

write.table(repr, "~/DATA/phylogeny/stepwise_taxon_repr.tsv", sep="\t", quote=FALSE, row.names = FALSE)  

         