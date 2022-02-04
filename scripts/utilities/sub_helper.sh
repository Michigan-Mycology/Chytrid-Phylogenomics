#!/bin/bash

LOGFILE=/scratch/tyjames_root/tyjames/amsesk/2021_neozygites/phylogeny/434_best_hits/iqtree_gene_trees/sub_helper.log
BATCHDIR=/scratch/tyjames_root/tyjames/amsesk/2021_neozygites/phylogeny/434_best_hits/iqtree_gene_trees/batch
MAX_CONCURRENT=65
IFS=$'\n'

function check_queue () {
    while :
    do
        queue=($(squeue -u $USER --format="%j" 2>&1))
        # Check for timeout error in squeue, restart loop if we find one 
        queue_timeout_error=$(echo ${queue[@]} | grep -i "socket" | wc -l)
        if [ $queue_timeout_error -gt 0 ]; then
            echo "Found error." > $LOGFILE
            sleep 60
            continue
        fi
        iq_in_queue=()
        for e in ${queue[@]}
        do
            case $e in
                gt_*)
                    iq_in_queue+=(${e})
                    ;;
            esac
        done
        echo ${#iq_in_queue[@]}
        break
    done
}

echo "IQTree Submission Helper Script Log" > $LOGFILE
echo "Maximum Concurrent Jobs: $MAX_CONCURRENT" >> $LOGFILE
echo "--------------------------------------" >> $LOGFILE
while read p; do
    
    #job_script=$(echo $p | cut -f1).sh
    marker=$(echo $p | cut -f2)
    job_script=${marker}.sh
    iq_in_queue=$(check_queue)

    echo $iq_in_queue
    while [ $iq_in_queue -eq  $MAX_CONCURRENT ]; do
        echo "$(date) ::: Too many jobs running right now. Checking again in 1 minute." >> $LOGFILE
        
        sleep 60
        
        iq_in_queue=$(check_queue)
    done    
    
    echo "$(date) ::: Room in queue. Submiting ${job_script} for ${marker}" >> $LOGFILE
    sbatch ${BATCHDIR}/${job_script}
    sleep 1

done < /scratch/tyjames_root/tyjames/amsesk/2021_neozygites/phylogeny/434_best_hits/iqtree_gene_trees/iqtree_job_mapper.tsv
