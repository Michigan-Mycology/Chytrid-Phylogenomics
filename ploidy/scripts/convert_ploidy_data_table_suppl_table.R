library(tidyverse)
library(readxl)

pdt = read_delim("~/DATA/pursuit/ploidy_final/ploidy_data_table.csv", delim = ",", col_names = F)
colnames(pdt) = c("Strain", "Assembly", "R1", "R2")

pdt = pdt %>%
  dplyr::select(-Assembly) %>%
  mutate(R1 = basename(R1)) %>%
  mutate(R2 = basename(R2))

# Remove strains sequenced in this paper
pdt = pdt[-seq(1,46,1),]

# Remove Coelomomyces, Rallo, Rmulti, PSC023
pdt = pdt[-c(8, 39, 40, 41, 42),]

# Remove fastq suffixes
pdt = pdt %>%
  mutate(R1 = gsub("[_]1[.]fastq", "", R1)) %>%
  mutate(R2 = gsub("[_]2[.]fastq", "", R1)) %>%
  dplyr::rename(`SRA Run ID (SRR)` = R1) %>%
  dplyr::select(-R2)

# Add SRA numbers for those missing
pdt[33,2] = "ERR2451319"
pdt[34,2] = "ERR2451316"
pdt[35,2] = "ERR2451317"
pdt[36,2] = "ERR2451318"

# Replace `_` with `-` in ranged SRR numbers
pdt = pdt %>%
  mutate(`SRA Run ID (SRR)` = gsub("[_]", "-", `SRA Run ID (SRR)`))

# Replace ploidy file names with full tree names
isolates = read_xlsx("~/work/pursuit/sheets/Pursuit_Isolates.xlsx", sheet = 1) %>%
  dplyr::select(SPECIES.TREE.LABEL, ploidy_file_prefix)

new_labels = c()
for (i in pdt$Strain) {
  new_labels = c(new_labels, isolates %>% filter(ploidy_file_prefix == i) %>% .$SPECIES.TREE.LABEL)
}

pdt$Strain = new_labels
pdt = pdt %>%
  mutate(Strain = gsub("[_]", " ", Strain))

View(pdt)

write.table(pdt, file = "~/work/pursuit/sheets/ploidy_srr_reference.tsv", sep = "\t", quote = F, row.names = F)
