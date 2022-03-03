### Suffixes and Such ###
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

		-s|--snpstats)
		SNPSTATS="$2"
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
echo "SNPSTATS: ${SNPSTATS}"

### Run the batch commands
while read p; do
	strain=`echo $p | cut -d ',' -f1`
	asm=`basename $(echo $p | cut -d ',' -f2)`
	python ~/work/Chytrid-Phylogenomics/ploidy/scripts/snp_stats_withPos_to_window_counts.py --strain ${strain} --l50fasta ${L50DIR}/${asm}.l50 --snp_stats ${SNPSTATS}/${strain}.snp_stats.tsv
done < $DATATABLE
