#!/bin/bash
#SBATCH --job-name=samtools
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=12g
#SBATCH --time=100:00:00
#SBATCH --account=tyjames1
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
        echo "Skippping ${f}"
        continue
    fi

    samtools sort --output-fmt BAM -o ${g} -@ 12 -m 10G ${f}
    #PicardCommandLine SortSam INPUT=$f OUTPUT=$g SORT_ORDER=coordinate
    
    samtools rmdup -s --output-fmt BAM ${g} ${h}
    #PicardCommandLine MarkDuplicates INPUT=$g OUTPUT=$h METRICS_FILE=$i

    samtools addreplacerg -r ID:S1 -r LB:L1 -r SM:${counter} --output-fmt BAM -@ 12 -o ${j} ${h}
    #PicardCommandLine AddOrReplaceReadGroups INPUT=$h OUTPUT=$j RGID=$counter RGLB=lib1 RGPU=1 RGPL=pacbio RGSM=$l && 
    
    samtools index ${j}
    #PicardCommandLine BuildBamIndex INPUT=$j
    
    counter=$((counter+1))
done
