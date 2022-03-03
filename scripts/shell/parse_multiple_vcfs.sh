#!/bin/zsh

L50_HEADERS_DIR=/home/aimzez/DATA/pursuit/ploidy/l50_headers
L50_HEADERS_EXT=l50_headers
KMERCOUNTS_DIR=/home/aimzez/DATA/pursuit/ploidy/kmerhist

while read p; do
	asm=`echo $p | cut -d ',' -f2 | xargs -L 1 -I @ basename @`
	strain=`echo $p | cut -d ',' -f1`
	
	if [ -f ${strain}.pdf ]; then
		echo "Skipping ${strain} because a PDF already exists."
		continue
	fi

	echo "Starting ${strain}"

	python ~/work/Chytrid-Phylogenomics/scripts/python/vcf_to_af.py \
		--strain ${strain} \
		-c ${L50_HEADERS_DIR}/${asm}.${L50_HEADERS_EXT} \
		--skip_na_mqrs \
		--dpfilt \
		--mqrsfilt \
		--plot \
		--kmercounts ${KMERCOUNTS_DIR}/${strain}_23 \
		${strain}.vcf

	echo "Done with ${strain}"

done < $1
