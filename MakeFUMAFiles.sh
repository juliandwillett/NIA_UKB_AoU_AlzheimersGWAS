#!/bin/bash

#SBATCH -n 1
#SBATCH -N 1
#SBATCH -J prep_fuma                                     # A single job name for the array
#SBATCH -t 0-01:00                                                   # Runtime in D-HH:MM
#SBATCH -p serial_requeue                                            # Partition to submit to
#SBATCH --mem=20000                                                  # Memory pool for all cores (see also --mem-per-cpu). 60G
#SBATCH -o extract_favor/query%A.out                                 # Standard output
#SBATCH -e extract_favor/query%A.err                                 # Standard error
#SBATCH --array=1-22                                                # Range of the array  
#SBATCH --open-mode=append                                           # Instructions to bypass parts of work that have already been completed.

echo ${SLURM_ARRAY_TASK_ID}

cd /n/home09/jwillett/true_lab_storage/00_AoU/
head -1 aou_ukb_allvar_meta_het_freq_qc_sorted_nomcc.txt > puma_files/meta_chr${SLURM_ARRAY_TASK_ID}.txt
awk -v var=$SLURM_ARRAY_TASK_ID '$13 == var {print $0}' aou_ukb_allvar_meta_het_freq_qc_sorted_nomcc.txt >> puma_files/meta_chr4.txt
gzip -9 puma_files/meta_chr${SLURM_ARRAY_TASK_ID}.txt
