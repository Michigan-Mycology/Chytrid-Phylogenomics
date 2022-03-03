library(tidyverse)

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  rename(assembly = `teamchytrid:Assembly/Final/genomes`) %>%
  mutate(assembly = gsub("[.]gz$", "", assembly)) %>%
  rename(proteome = `teamchytrid:Phylogeny/Sequences/pep`) %>%
  mutate(proteome = gsub("[.]gz$", "", proteome))

asmstats = read_delim("~/work/Chytrid-Phylogenomics/figures/Supplmental/all_assembly_stats.tsv", delim='\t', col_names = F) %>%
  rename(assembly = X1, asmlen = X2, gc = X3, n50 = X4, l50 = X5) %>%
  full_join(isolates %>% select(SPECIES.TREE.LABEL, assembly)) %>%
  filter(!is.na(SPECIES.TREE.LABEL)) %>%
  filter(!SPECIES.TREE.LABEL %in% c("Rozella_rhizoclosmatii",
                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Amber.LCG",
                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Orange.LCG"))

nprots = read_delim("~/work/Chytrid-Phylogenomics/figures/Supplmental/all_pep_protein_numbers.tsv", delim="\t", col_names = F) %>%
  rename(proteome = X1, nprots = X2)  %>%
  full_join(isolates %>% select(SPECIES.TREE.LABEL, proteome, `this project`)) %>%
  filter(!SPECIES.TREE.LABEL %in% c("Rozella_rhizoclosmatii",
                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Amber.LCG",
                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Orange.LCG"))

busco_protein = read_delim("~/work/Chytrid-Phylogenomics/figures/Supplmental/all_busco_prot.tsv", delim=",", col_names = F) %>%
  rename(proteome = X1, complete_buscos = X2, complete_and_duplicated_buscos = X3) %>%
  mutate(total_buscos_in_db = 758) %>%
  mutate(busco_completeness = complete_buscos/total_buscos_in_db) %>%
  mutate(busco_duplication = complete_and_duplicated_buscos/complete_buscos) %>%
  full_join(isolates %>% select(SPECIES.TREE.LABEL, proteome)) %>%
  filter(!is.na(SPECIES.TREE.LABEL)) %>%
  filter(!SPECIES.TREE.LABEL %in% c("Rozella_rhizoclosmatii",
                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Amber.LCG",
                                    "Coelomomyces_lativittatus_CIRM-AVA-1-Orange.LCG"))

combined_table = asmstats %>%
  left_join(nprots) %>%
  left_join(busco_protein %>% select(SPECIES.TREE.LABEL, busco_completeness, busco_duplication)) %>%
  select(-assembly, -proteome) %>%
  select(SPECIES.TREE.LABEL, `this project`, asmlen, gc, n50, l50, nprots, busco_completeness, busco_duplication) %>%
  arrange(SPECIES.TREE.LABEL)

write.table(combined_table,
            file = "~/work/Chytrid-Phylogenomics/figures/Supplmental/Supplementary_Table_1_genomes_summary.tsv",
            quote = F,
            row.names = F,
            sep="\t")
