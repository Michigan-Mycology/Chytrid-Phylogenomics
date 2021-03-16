#!/bin/bash

L50SUFF="l50"

### Parse Command Line Arguments ###
POSITIONAL=()
while [[ $# -gt 0 ]]; do

	key="$1"

	case $key in
		
		-l|--l50dir)
		L50DIR="$2"
		shift # past argument
		shift # past value
		;;

		-d|--datatable)
		DATATABLE="$2"
		shift
		shift
		;;

		*)    # unknown option
    		POSITIONAL+=("$1") # save it in an array for later
    		shift # past argument
		;;
	esac
done

set -- "${POSITIONAL[@]}" # restore positional parameters

### Print Argument Values before Running
echo "L50DIR: ${L50DIR}"
echo "DATATABLE: ${DATATABLE}"

while read p; do
    strain=$(echo $p | cut -d ',' -f1)
    l50_asm=$(basename `echo $p | cut -d ',' -f2`).${L50SUFF}
    snp_stats=${strain}.snp_stats.tsv
    out=${strain}.counts.tsv
    python ~/dev/Chytrid-Phylogenomics/ploidy/scripts/snp_contig_counts.py -a ${L50DIR}/${l50_asm} -o ${out} ${snp_stats}
done < $DATATABLE
