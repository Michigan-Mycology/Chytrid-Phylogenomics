library(tidyverse)
library(readxl)
library(harrypotter)

mat = read_delim("~/DATA/bonitomes/topology_heatmatp.tsv", delim="\t")
mat = mat %>%
  select(-grep("[(][0-9], [(][0-9], [0-9][)][)]", colnames(mat)))
mat = mat %>%
  pivot_longer(cols = colnames(mat)[-1], names_to = "sister", values_to = "count") %>%
  rename(group = `...1`) %>%
  mutate(count_log = log(count+1))
mat = mat %>%
  mutate(sister = factor(sister, levels = c("0", "1", "2", "3", "(0, 1)", "(0, 2)", "(0, 3)", "(1, 2)", "(1, 3)", "(2, 3)"))) %>%
  mutate(group = factor(group, levels = c("0", "1", "2", "3", "(0, 1)", "(0, 2)", "(0, 3)", "(1, 2)", "(1, 3)", "(2, 3)", "(0, (1, 2))", "(0, (1, 3))", "(0, (2, 3))", "(1, (0, 2))", "(1, (0, 3))", "(1, (2, 3))",
                                          "(2, (0, 1))", "(2, (0, 3))", "(2, (1, 3))", "(3, (0, 1))", "(3, (0, 2))", "(3, (1, 2))")))

hm = ggplot(mat, aes(x = group, y=sister, fill = count_log)) +
  geom_tile(color="white") +
  geom_text(aes(label=count)) +
  scale_fill_hp(house = "ravenclaw", na.value = "grey88") +
  xlab("Sister 1") +
  ylab("Sister 2") +
  coord_fixed() +
  theme_bw() +
  theme(
    axis.text.x = element_text(angle = 90, hjust =1)
  )
hm

ggsave(filename = "~/dev/EHB-Bonitomes/Intermediate_Figures/trees/gt_topology_counts_heatmap.pdf",
       plot = hm,
       width = 8.5,
       height = 11,
       device="pdf",
       unit = "in")
ggsave(filename = "~/dev/EHB-Bonitomes/Intermediate_Figures/trees/gt_topology_counts_heatmap.png",
       plot = hm,
       width = 8.5,
       height = 11,
       device="png",
       unit = "in")
