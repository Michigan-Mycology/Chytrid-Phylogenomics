library(phytools)
library(tidyverse)
library(readxl)
library(magrittr)
library(scales)

##########################
### Mono/Poly/Missing ####
##########################

args = commandArgs(trailingOnly = TRUE)
matrix_path = args[1]
isolates_xlsx = args[2]

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

named_mat$value = factor(named_mat$value, levels=c("-1","0","1"))
named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

replaced = named_mat
replaced$value = recode(named_mat$value, "Missing" = 0, "Polyphyletic" = 1, "Monophyletic" = 1)
replaced %<>%
  pivot_wider(names_from = "marker", values_from="value")

isolate_repr = replaced %>%
  mutate(rowsum = rowSums(.[2:dim(replaced)[2]])) %>%
  mutate(rowsum = rowsum/dim(replaced)[2]) %>%
  select(SPECIES.TREE.LABEL, rowsum) %>%
  rename(Trimal_Phyly_Score_ManualCuration_SpikeIn = rowsum) %>%
  arrange(SPECIES.TREE.LABEL)
write.table(isolate_repr, "isolate_repr.tsv", sep="\t", row.names = FALSE, quote=FALSE)
write.table(replaced, "occupancy_matrix.tsv", sep="\t", row.names = FALSE, quote=FALSE)

colors = c("grey", muted("red"), muted("blue"))
names(colors) = levels(named_mat$value)
hm = ggplot(named_mat, aes(x=marker, y=SPECIES.TREE.LABEL, fill=value)) + 
  geom_tile() +
  scale_fill_manual(values=colors) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=4)
  )
ggsave("phyly_heatmap.pdf", hm, width = unit(17, "in"), height = unit(22, "in"))