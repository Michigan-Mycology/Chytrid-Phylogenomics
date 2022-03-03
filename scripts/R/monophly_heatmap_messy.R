library(tidyverse)
library(readxl)
library(scales)
library(magrittr)

##########################
### Mono/Poly/Missing ####
##########################

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
#p_monophyletic <- function(column) {
#  count = length(column[(column == 1 | column == 0)])
#  count / length(column)
#}

#replaced = named_mat
#replaced$value = recode(named_mat$value, "Missing" = 0, "Polyphyletic" = 1, "Monophyletic" = 1)
#replaced %<>%
#  pivot_wider(names_from = "marker", values_from="value")

#marker_occ = replaced %>%
#  select(-SPECIES.TREE.LABEL) %>% 
#  summarise_all("sum") %>%
#  mutate_all(function (x) x/140) %>%
#  gather
#summary(marker_occ$value)
#hist(marker_occ$value)

'
isolate_repr = replaced %>%
  mutate(rowsum = rowSums(.[2:500])) %>%
  mutate(rowsum = rowsum/499) %>%
  select(SPECIES.TREE.LABEL, rowsum) %>%
  rename(Round4_PostSpike = rowsum) %>%
  arrange(SPECIES.TREE.LABEL)
write.table(isolate_repr, "~/DATA/phylogeny_trimfix/round4_mancur_round3/isolate_repr.tsv", sep="\t", row.names = FALSE, quote=FALSE)


stepwise_taxon_repr = isolate_repr
stepwise_taxon_repr = stepwise_taxon_repr %>%
  left_join(isolate_repr, by = "SPECIES.TREE.LABEL")
stepwise_taxon_repr = stepwise_taxon_repr %>%
  left_join(isolate_repr, by = "SPECIES.TREE.LABEL")
raw_repr = read_delim("~/DATA/phylogeny/taxon_repr_pre_trimal.tsv", delim="\t") %>%
  mutate(Raw = occurrences_in_unaln / 582) %>%
  select(SPECIES.TREE.LABEL, Raw)
stepwise_taxon_repr = stepwise_taxon_repr %>%
  left_join(raw_repr, by = "SPECIES.TREE.LABEL")

stepwise_taxon_repr %<>% 
  select(SPECIES.TREE.LABEL,Raw, Trimal, Trimal_Phyly, Trimal_Phyly_Score) %>%
  arrange(SPECIES.TREE.LABEL)
write.table(stepwise_taxon_repr, "stepwise_taxon_repr.tsv", sep="\t", col.names=TRUE, row.names = FALSE, quote=FALSE)

## Marker occupancy now
count_missing <- function(column) {
  count = 0
  for (r in column) {
    if (r == -1) {
      count = count + 1
    }
  }
  1- (count/length(column))
}
hist = mat %>% summarise_all(count_missing)
hist(as.numeric(hist[1,]), breaks=25, xlim=c(0.6,1.0), main="Marker Occupancy", xlab="Occupancy")
abline(v=0.75, col="red")

###############
### nLeaves ###
###############
mat = read_delim("~/DATA/phylogeny/leaves_of_common_ancestor.csv", delim=",")
colnames(mat)[1] = "LTP"
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP)
named_mat = mat %>% 
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP) %>%
  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%

#ord = hclust(dist(named_mat[,-1], method="euclidean"), method="ward.D")$ord


#named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

ggplot(named_mat, aes(x=marker, y=SPECIES.TREE.LABEL, fill=value)) + 
  geom_tile() +
  scale_fill_continuous(type="viridis") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=4)
  )

#########################
### duplication level ###
#########################
mat = read_delim("~/DATA/phylogeny/duplication_level.csv", delim=",")
colnames(mat)[1] = "LTP"
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP)
named_mat = mat %>% 
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP) %>%
  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%
  mutate(value = log(value))

#ord = hclust(dist(named_mat[,-1], method="euclidean"), method="ward.D")$ord


#named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

ggplot(named_mat, aes(x=marker, y=SPECIES.TREE.LABEL, fill=value)) + 
  geom_tile() +
  scale_fill_continuous(type="viridis") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=4)
  )

#########################
### max score difference ###
#########################
mat = read_delim("~/DATA/phylogeny/max_score_difference.csv", delim=",")
colnames(mat)[1] = "LTP"
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, LTP)
named_mat = mat %>% 
  full_join(isolates) %>%
  select(SPECIES.TREE.LABEL, LTP, everything()) %>%
  select(-LTP) %>%
  #select(-high_mean) %>%
  gather(key=marker, "value", -SPECIES.TREE.LABEL) %>%
  mutate(value = sqrt(value))

#ord = hclust(dist(named_mat[,-1], method="euclidean"), method="ward.D")$ord


#named_mat$value = recode(named_mat$value, "-1" = "Missing", "0" = "Polyphyletic", "1" = "Monophyletic")

ggplot(named_mat, aes(x=marker, y=SPECIES.TREE.LABEL, fill=value)) + 
  geom_tile() +
  #scale_fill_continuous(type="viridis")+
  scale_fill_gradient(low = "royalblue1", high = muted("red"), na.value="white") +
  theme(
    axis.text.x = element_text(angle=45, hjust=1, size=4)
  )

high_mean = mat %>%
  select(-LTP) %>%
  mutate_all(sqrt) %>%
  summarise_all(mean, na.rm = TRUE) %>%
  gather %>%
  filter (value >= 10) %>%
  .$key
named_mat %>% select(-high_mean)
###################
d = read_delim("~/DATA/phylogeny/hit_report_all.csv", delim="\t", col_names = FALSE) %>%
  mutate(X3 = as.numeric(X3))

plot(log(d$X3), d$X4, cex=0.2)
'