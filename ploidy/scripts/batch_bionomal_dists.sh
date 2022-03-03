#!/bin/zsh

CHYTRID_PHYLO=/home/aimzez/work/Chytrid-Phylogenomics
DEPTHS=/home/aimzez/DATA/pursuit/ploidy_final/contig.depth.l50
SNP_STATS=/home/aimzez/DATA/pursuit/ploidy_final/snp_stats

for i in $(ls $DEPTHS | sed 's/.sorted.bam.depth.ContigMean.L50//'); do
	Rscript ${CHYTRID_PHYLO}/ploidy/scripts/bionomial_dists.R ${DEPTHS}/${i}.sorted.bam.depth.ContigMean.L50 ${SNP_STATS}/${i}.snp_stats.tsv
done
