#!/bin/bash
# step1_2
#SBATCH -J step1_2         # Job name: step1_2
#SBATCH -p shared          # Use the shared partition
#SBATCH -c 4               # Number of CPU cores: 4
#SBATCH -t 3-00:00         # Time limit: 3 days
#SBATCH --mem 4000         # Memory limit: 4000 MB
#SBATCH -o step1_2.out     # File for standard output
#SBATCH -e step1_2.err     # File for standard error

# Defining directories and files on DNAnexus
pfile_dir="project-GKGxZGjJ15P9gg1J6Z7q4ZzF:/MW/WGS_chunks/WGS_runs/biallelic_pfiles" # Directory for PGEN files
pheno_dir="project-GKGxZGjJ15P9gg1J6Z7q4ZzF:/MW/Pheno"                              # Directory for phenotype files
pred_dir="project-GKGxZGjJ15P9gg1J6Z7q4ZzF:/MW/Step1"                               # Directory for prediction files
pred="ukb_step1_BT_pred.list"                                                      # Prediction file list
pheno_file="pheno2"                                                                # Phenotype file
covar_file="covar2"                                                                # Covariate file
wgs_reg="project-GKGxZGjJ15P9gg1J6Z7q4ZzF:/MW/WGS_chunks/WGS_runs/Step1_2_corrected" # Output directory

list="/n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/pgens.txt"                         # List of PGEN prefixes

# Processing each line in the list of PGEN files
while read -r line; do
    prefix="$line"

    # Command to run regenie for step 2 analysis
    reg_cmd="regenie \
    --step 2 \
    --pgen ${prefix} \
    --phenoFile ${pheno_file} \
    --covarFile ${covar_file} \
    --catCovarList Shipment_batch \
    --maxCatLevels 25 \
    --bt \
    --firth --firth-se --approx --pThresh 0.05 \
    --minMAC 3 \
    --lowmem \
    --threads 2 \
    --pred ${pred} \
    --bsize 400 \
    --out 1+2_${prefix}"

    # Running the command on DNAnexus using Swiss Army Knife app
    dx run app-swiss-army-knife \
    -iin="${pfile_dir}/${prefix}.pgen" \
    -iin="${pfile_dir}/${prefix}.pvar" \
    -iin="${pfile_dir}/${prefix}.psam" \
    -iin="${pred_dir}/${pred}" \
    -iin="${pred_dir}/ukb_step1_BT_1.loco" \
    -iin="${pred_dir}/ukb_step1_BT_2.loco" \
    -iin="${pheno_dir}/${pheno_file}" \
    -iin="${pheno_dir}/${covar_file}" \
    -icmd="${reg_cmd}" \
    --tag="WGS_step1_2" \
    --instance-type "mem1_ssd1_v2_x2" \
    --name "WGS_step1_2_chunk" \
    --destination="${wgs_reg}" \
    --yes

done < "${list}"
