library(tidyverse)

clust1 = read_delim("~/data/candida/snpeff/cluster_snps/corrected/cluster1_snps_CORRECTED.csv", delim = "\t") %>%
  select(-"...1") %>%
  rename("833-4" = "14", "838-4" = "45", "833-21" = "3", "833-27" = "44", "833-19" = "8", "838-2" = "39") %>%
  #filter_at(vars(`833-27`, `833-19`, `838-2`, `833-4`, `838-4`, `833-21`), all_vars(. != "N")) %>%
  filter( (`833-19` == `833-21`) & (`833-19` != `833-4` & `833-19` !=`838-4` & `833-19` != `833-27` & `833-19` !=`838-2`) )

write_delim(clust1,
            file = "~/data/candida/snpeff/cluster_snps/corrected/cluster1_833-19.833-21_shared_unique_snps_CORRECTED.tsv",
            delim = "\t",
            quote_escape = "none")

clust2 = read_delim("~/data/candida/snpeff/cluster_snps/corrected/cluster2_snps_CORRECTED.csv", delim = "\t") %>%
  select(-"...1") %>%
  rename("814-42" = "11", "814-147" = "12", "814-168" = "19", "814-45" = "17", "814-183" = "25", "814-186" = "37") %>%
  filter_at(vars(`814-42`, `814-147`, `814-168`, `814-45`, `814-183`, `814-186`), all_vars(. != "N")) %>%
  filter(`814-168` != `814-42` & `814-168` !=`814-147` & `814-168` != `814-45` & `814-168` !=`814-183` & `814-168` != `814-186`)

write_delim(clust2,
            file = "~/data/candida/snpeff/cluster_snps/cluster2_814-168_unique_snps_CORRECTED.tsv",
            delim = "\t",
            quote_escape = "none")


clust2_45_42 = read_delim("~/data/candida/snpeff/cluster_snps/corrected/cluster2_snps_CORRECTED.csv", delim = "\t") %>%
  select(-"...1") %>%
  rename("814-42" = "11", "814-147" = "12", "814-168" = "19", "814-45" = "17", "814-183" = "25", "814-186" = "37") %>%
  filter_at(vars(`814-42`, `814-147`, `814-168`, `814-45`, `814-183`, `814-186`), all_vars(. != "N")) %>%
  filter(`814-42` == `814-45` &
           `814-42` !=`814-168` &
           `814-42` != `814-147` &
           `814-42` !=`814-183` &
           `814-42` != `814-186`)

write_delim(clust2,
            file = "~/data/candida/snpeff/cluster_snps/cluster2_814-168_unique_snps_CORRECTED.tsv",
            delim = "\t",
            quote_escape = "none")
