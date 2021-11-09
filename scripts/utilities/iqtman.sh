#!/bin/bash
#SBATCH --job-name=iqtman
#SBATCH --mail-type=BEGIN,END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=250m
#SBATCH --time=336:00:00
#SBATCH --account=tyjames1
#SBATCH --partition=standard

sh /scratch/tyjames_root/tyjames/amsesk/2021_neozygites/phylogeny/434_best_hits/iqtree_gene_trees/sub_helper.sh
