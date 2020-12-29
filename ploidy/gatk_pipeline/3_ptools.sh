#!/bin/bash
#SBATCH --job-name=ptools
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=25g
#SBATCH --time=100:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard

module load Bioinformatics
module load samtools
module load picard-tools

FFILES=$1

counter=1
for f in $FFILES/*.sam
do
    g=${f/.sam/.sorted.bam}
    h=${f/.sam/.dedupped.sorted.bam}
    i=${f/.sam/.duplicate.metrics.txt}
    j=${f/.sam/.dedupped.sorted.addrg.bam}
    k=$(basename $f)
    l=${k/.sam/}
  
    if ls $j; then
        counter=$((counter+1))
        continue
    fi

    PicardCommandLine SortSam INPUT=$f OUTPUT=$g SORT_ORDER=coordinate
    PicardCommandLine MarkDuplicates INPUT=$g OUTPUT=$h METRICS_FILE=$i && \
    PicardCommandLine AddOrReplaceReadGroups INPUT=$h OUTPUT=$j RGID=$counter RGLB=lib1 RGPU=1 RGPL=pacbio RGSM=$l && \
    PicardCommandLine BuildBamIndex INPUT=$j
    
    counter=$((counter+1))
done
