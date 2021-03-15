#!/bin/bash

L50DIR=/scratch/aimzez/pursuit/data/assemblies_l50


while read p; do
    strain=$(echo $p | cut -d ',' -f1)
    l50_asm=$(basename `echo $p | cut -d ',' -f2`).l50
    snp_stats=${strain}.snp_stats.tsv
    out=${strain}.counts.tsv
    python ~/dev/Chytrid-Phylogenomics/ploidy/scripts/snp_contig_counts.py -a ${L50DIR}/${l50_asm} -o ${out} ${snp_stats}
done < $1
