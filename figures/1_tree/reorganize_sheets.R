library(readxl)
library(tidyverse)

#### Merge ploidy file prefixes
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, `teamchytrid:Assembly/Final/genomes`) %>%
  rename(assembly = `teamchytrid:Assembly/Final/genomes`) %>%
  mutate(assembly = gsub("[.]gz", "", assembly))

orig_pdt = read_delim("~/DATA/pursuit/ploidy_final/ploidy_data_table.csv", delim = ',', col_names = F)
pdt = read_delim("~/DATA/pursuit/ploidy_final/ploidy_data_table.csv", delim = ',', col_names = F) %>%
  rename(ploidy_file_prefix = X1, assembly = X2, R1 = X3, R2 = X4) %>%
  mutate(assembly = basename(assembly)) %>%
  select(-R1, -R2) %>%
  left_join(isolates, by="assembly")

isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(-ploidy_file_prefix) %>%
  left_join(pdt %>% select(SPECIES.TREE.LABEL, ploidy_file_prefix))

write.table(isolates, file= "~/work/pursuit/sheets/isolates-tmp.tsv", sep="\t", quote=F, row.names = F)

#### Merge rabern's gene stuff into growing trait table
ultra = read_delim("~/work/pursuit/sheets/Phylogenomics_Ultrastructure_renamed.tsv", delim="\t")
more_rabern = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates_Phylo_Traits.xlsx") %>%
  select(SPECIES.TREE.LABEL, EFL:Selenocystein) %>%
  mutate_at(vars(-SPECIES.TREE.LABEL), ~ifelse(is.na(.), 0, 1))


merged = ultra %>%
  left_join(more_rabern)

write.table(x=merged, file="~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv",
            sep="\t", row.names = F, col.names = T, quote = F)


#### Merge Final Ploidy Calls into growing trait table

ploidy = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, coding, `final call of life cycle`)

next_merged = merged %>%
  left_join(ploidy)
next_merged

write.table(x=next_merged, file="~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv",
            sep="\t", row.names = F, col.names = T, quote = F)

#### Add L50 assembly length to growing Phylo Traits table
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx") %>%
  select(SPECIES.TREE.LABEL, `teamchytrid:Assembly/Final/genomes`) %>%
  rename(assembly = `teamchytrid:Assembly/Final/genomes`) %>%
  mutate(assembly = gsub("[.]gz", "", assembly))

phyltraits = read_delim('~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv', delim="\t")

l50_asmlen = read_delim("~/DATA/pursuit/ploidy_final/all.l50.asmlen", delim="\t", col_names = F) %>%
  rename(assembly = X1, l50_assembly_length = X2) %>%
  mutate(assembly = gsub("[.]l50", "", assembly)) %>%
  left_join(isolates)

phyltraits_new = phyltraits %>%
  left_join(l50_asmlen)

write.table(x=phyltraits_new, file="~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv",
            sep="\t", row.names = F, col.names = T, quote = F)


#### Did we do ploidy? Add that growing table
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", na="NA") %>%
  select(SPECIES.TREE.LABEL, ploidy_file_prefix)

phyltraits = read_delim('~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv', delim="\t", na = "NA")

phyltraits = phyltraits %>%
  left_join(isolates) %>%
  rename(UM_ploidy = ploidy_file_prefix) %>%
  mutate(UM_ploidy = ifelse(is.na(UM_ploidy), 0, 1))

write.table(x=phyltraits, file="~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv",
            sep="\t", row.names = F, col.names = T, quote = F)

#### Add cell cycle genes to growing table
cell_cycle = read_xlsx("~/work/pursuit/sheets/chytrid-G1Sregulators-4column.xlsx") %>%
  filter(!`Genus Species` %in% c("Non-fungal eukaryota (outgroup)", "Mucormycota ingroup", "Dikarya ingroup")) %>%
  separate(`Genus Species`, into = c("genus", "species", "strain"), extra="merge", sep=" ") %>%
  mutate(strain = stringr::str_replace(strain, "[(]", "")) %>%
  mutate(strain = stringr::str_replace(strain, "[)]", "")) %>%
  unite("species_name", genus, species, sep = " ") %>%
  mutate(species_name = stringr::str_trim(species_name, side="both"))

phyltraits = read_delim("~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv", delim="\t") %>%
  mutate(label_merge_with_nicks = SPECIES.TREE.LABEL) %>%
  separate(label_merge_with_nicks, into = c('genus', 'species', 'strain'), extra = 'merge', sep = "_") %>%
  unite("species_name", genus, species, sep=" ") %>%
  mutate(strain = stringr::str_replace(strain, "[.]*LCG", "")) %>%
  mutate(strain = stringr::str_replace(strain, "[_.]*v[0-9][.]*[0-9]*", "")) %>%
  mutate(strain = ifelse(strain %in% cell_cycle$strain, strain, NA))


combined = phyltraits %>%
  left_join(cell_cycle, by=c("species_name", "strain"))

View(
  combined %>%
    filter_at(vars(Ancestral_combined:Whi5_Fungal), all_vars(is.na(.))) %>%
    select(SPECIES.TREE.LABEL, species_name, strain, Ancestral_combined:Whi5_Fungal)
)

combined %>%
  filter_at(vars(Ancestral_combined:Whi5_Fungal), all_vars(is.na(.)) )

combined = combined %>%
  select(-species_name, -strain)

write.table(x=combined, file="~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv",
            sep="\t", row.names = F, col.names = T, quote = F)

#### merge in Rabern's new ef1a, etc - 041221
updated = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates_Phylo_Traits_Rs_041221.xlsx") %>%
  select(SPECIES.TREE.LABEL, EFL:RNR) %>%
  mutate_at(vars(EFL:RNR), ~ifelse(is.na(.), 0, 1))
original = read_delim("~/work/pursuit/sheets/Pursuit_Phylo_Traits_032621.tsv", delim="\t") %>%
  select(-EFL, -EF1a, -CblA, -CblC, -CblD, -MeaB, -MetH, -MMCoAepi, -MMCoAmut, -RNR)

full_updated = original %>%
  left_join(updated)

write.table(full_updated,
            file = "~/work/pursuit/sheets/Pursuit_Phylo_Traits_041221.tsv",
            quote = F,
            row.names = F,
            sep = "\t")
