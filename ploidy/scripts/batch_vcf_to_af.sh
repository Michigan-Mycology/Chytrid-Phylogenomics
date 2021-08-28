### Suffixes and Such ###
VCFSUFF="g.vcf"
L50HEADERSUFF="l50_headers"
CHYTRIDPHYLO=/home/aimzez/work/Chytrid-Phylogenomics/

### Parse Command Line Arguments ###
POSITIONAL=()
while [[ $# -gt 0 ]]; do

	key="$1"

	case $key in
		-k|--kmerhist)
    		KMERHIST="$2"
    		shift # past argument
    		shift # past value
    		;;
    		
		-l|--l50headers)
		L50HEADERS="$2"
		shift # past argument
		shift # past value
		;;

		-g|--gvcfs)
		GVCFS="$2"
		shift
		shift
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
echo "KMERHIST: ${KMERHIST}"
echo "L50HEADERS: ${L50HEADERS}"
echo "GVCFS: ${GVCFS}"
echo "DATATABLE: ${DATATABLE}"

### Run the batch commands
while read p; do
	asm=$(basename `echo $p | cut -d ',' -f2`)
	strain=`echo $p | cut -d ',' -f1`
	python ${CHYTRIDPHYLO}/ploidy/scripts/vcf_to_af.py --strain ${strain} -c ${L50HEADERS}/${asm}.${L50HEADERSUFF} --skip_na_mqrs --dpfilt --mqrsfilt --plot --kmercounts ${KMERHIST}/${strain}_23 ${GVCFS}/${strain}.${VCFSUFF}
done < $DATATABLE
