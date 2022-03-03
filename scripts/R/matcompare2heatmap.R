library(phytools)
library(tidyverse)
library(ggtree)
library(readxl)
library(magrittr)
library(RColorBrewer)
library(scales)

##########################
### Mono/Poly/Missing ####
##########################

args = commandArgs(trailingOnly = TRUE)
matrix_path = args[1]
isolates_xlsx = args[2]

#matrix_path = "~/DATA/phylogeny_trimfix/round1_raw/compare_matrix.csv"
#isolates_xlsx = "~/work/pursuit/sheets/Pursuit_Isolates.xlsx"

mat = read_delim(matrix_path, delim=",")

colnames(mat)[1] = "LTP"
isolates = read_xlsx(isolates_xlsx) %>%
  select(SPECIES.TREE.LABEL, LTP)
named_mat = mat %>% 
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP)

named_mat %<>%  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%
  mutate(value = as.character(value)) %>%
  mutate(value = as.factor(value))

named_mat$value = as.factor(named_mat$value)
named_mat$value = recode(named_mat$value, "11" = "Monophyletic -> Monophyletic",
                         "10" = "Monophyletic -> Polyphyletic",
                         "01" = "Polyphyletic -> Monophyletic",
                         "00" = "Polyphyletic -> Polyphyletic",
                         "-10" = "Missing -> Polyphyletic",
                         "-11" = "Missing -> Monophyletic",
                         "-1-1" = "Missing -> Missing",
                         "0-1" = "Polphyletic -> Missing",
                         "1-1" = "Monophyletic -> Missing"
                         )

cell_counts = as.data.frame(table(named_mat %>% select(-SPECIES.TREE.LABEL) %>% .$value))
write.table(x=cell_counts, file="cell_counts.tsv", quote=FALSE, row.names = FALSE, sep="\t")

hm = ggplot(named_mat, aes(x=marker, y=SPECIES.TREE.LABEL, fill=value)) + 
  geom_tile() +
  scale_fill_manual(values=c("grey", "black", muted("red"), "blue", "black", "red", muted("blue"))) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=4)
  )
ggsave("phyly_compare_heatmap.pdf", hm, width = unit(22, "in"), height = unit(17, "in"))
