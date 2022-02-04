#!/bin/bash
#SBATCH --job-name=samtools
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=15g
#SBATCH --time=168:00:00
#SBATCH --account=tromeara99
#SBATCH --partition=standard

module load Bioinformatics
module load samtools

FFILES=$1

counter=1
for f in $FFILES/*.sam
do
    g=${f/.sam/.sorted.bam}
    h=${f/.sam/.dedupped.sorted.bam}
    #i=${f/.sam/.duplicate.metrics.txt}
    j=${f/.sam/.dedupped.sorted.addrg.bam}
    k=$(basename $f)
    l=${k/.sam/}
  
    if [ -f $j ]; then
        echo Skippping ${f}$'\t'${counter}
        counter=$((counter+1))
        continue
    fi

    #samtools sort --output-fmt BAM -o ${g} -@ 12 -m 178G ${f}
    
    #samtools rmdup -s --output-fmt BAM ${g} ${h}

    #samtools addreplacerg -r ID:S1 -r LB:L1 -r SM:${counter} --output-fmt BAM -@ 12 -o ${j} ${h}
    
    #samtools index ${j}
    
    echo ${f}$'\t'${counter}
    counter=$((counter+1))
done
