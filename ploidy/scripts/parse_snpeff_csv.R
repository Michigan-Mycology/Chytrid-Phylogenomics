library(tidyverse)

clust1 = read_delim("~/DATA/candida/snpeff_cluster_snps/cluster1_snps.csv", delim = "\t") %>%
  select(-X1) %>%
  rename("833-4" = "11", "838-4" = "12", "833-21" = "13", "833-27" = "3", "833-19" = "5", "838-2" = "8") %>%
  #filter_at(vars(`833-27`, `833-19`, `838-2`, `833-4`, `838-4`, `833-21`), all_vars(. != "N")) %>%
  filter( (`833-19` == `833-21`) & (`833-19` != `833-4` & `833-19` !=`838-4` & `833-19` != `833-27` & `833-19` !=`838-2`) )

write_delim(clust1,
            file = "~/DATA/candida/snpeff_cluster_snps/cluster1_833-19.833-21_shared_unique_snps.tsv",
            delim = "\t", 
            quote_escape = "none")

clust2 = read_delim("~/DATA/candida/snpeff_cluster_snps/cluster2_snps.csv", delim = "\t") %>%
  select(-X1) %>% 
  rename("814-42" = "14", "814-147" = "17", "814-168" = "19", "814-45" = "23", "814-183" = "25", "814-186" = "37") %>%
  filter_at(vars(`814-42`, `814-147`, `814-168`, `814-45`, `814-183`, `814-186`), all_vars(. != "N")) %>%
  filter(`814-168` != `814-42` & `814-168` !=`814-147` & `814-168` != `814-45` & `814-168` !=`814-183` & `814-168` != `814-186`)

write_delim(clust2,
            file = "~/DATA/candida/snpeff_cluster_snps/cluster2_814-168_unique_snps.tsv",
            delim = "\t", 
            quote_escape = "none")
