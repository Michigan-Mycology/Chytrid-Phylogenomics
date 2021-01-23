#!/bin/bash
#SBATCH --job-name=asm_prep
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8g
#SBATCH --time=4:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard

module load Bioinformatics samtools picard-tools

ASMDIR=$1
for f in $(ls $ASMDIR/*.fasta)
do
    bwa index ${f}
    samtools faidx ${f}
    PicardCommandLine CreateSequenceDictionary R=${f}
done
