#!/bin/zsh

ASMDIR=/home/aimzez/DATA/pursuit/data/
VCFDIR=.

TABLE=$1
OUT=$2

echo "strain,assembly,assembly_length,snps,snp_density" > $OUT
while read p; do
	IFS=,
	spl=($p)
	strain=${spl[0]}
	assembly=`basename ${spl[1]}`
	parsed_vcf=${strain}.snp_stats.tsv
	python ~/work/Chytrid-Phylogenomics/scripts/ploidy/snp_density/calc_snp_density.py \
		--assembly ${ASMDIR}/${assembly} \
		--parsed_vcf ${VCFDIR}/${parsed_vcf} \
		--strain ${strain}
done < $TABLE >> $OUT
