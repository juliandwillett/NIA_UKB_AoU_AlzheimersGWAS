#script to generate Step 1 outputs

#!/bin/bash
# step1
#SBATCH -J step1
#SBATCH -p test
#SBATCH -c 24
#SBATCH -t 0-08:00
#SBATCH --mem 80000
#SBATCH -o pgen_%A_%a.out # Standard output
#SBATCH -e pgen_%A_%a.err # Standard error

regenie \
  --step 1 \
  --pgen /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/pre-Step1/ukb_all_chrs-merge \
  --extract /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/pre-Step1/qc_pass.snplist \
  --keep /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/pre-Step1/qc_pass.id \
  --phenoFile /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/pheno_files/pheno2 \
  --covarFile /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/pheno_files/covar2 \
  --catCovarList Shipment_batch \
  --maxCatLevels 25 \
  --bt \
  --bsize 1000 \
  --threads 24 \
  --out ukb_step1_BT
