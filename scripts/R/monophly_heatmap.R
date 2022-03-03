library(tidyverse)
library(readxl)
library(scales)
library(magrittr)

args = commandArgs(trailingOnly = TRUE)
mat = args[1]
ltp_mapper = args[2]

mat = read_delim(mat, delim=",")

colnames(mat)[1] = "LTP"
isolates = read_delim(ltp_mapper, delim='\t', col_names=FALSE) %>%
  rename(SPECIES.TREE.LABEL = X1, LTP = X2)
named_mat = mat %>% 
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP)

named_mat %<>%  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%
  mutate(value = as.character(value)) %>%
  mutate(value = as.factor(value))

named_mat$value = factor(named_mat$value, levels=c("-1","0","1"))
named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

hm = ggplot(named_mat, aes(x=marker, y=SPECIES.TREE.LABEL, fill=value)) + 
  geom_tile() +
  scale_fill_manual(values=c("grey", muted("red"), muted("blue"))) +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=4)
  )

ggsave(filename = "monophyly_heatmap.pdf", plot=hm, width=17, height=22)
