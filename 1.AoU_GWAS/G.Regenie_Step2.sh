# Get files for current chromosome
gsutil -m cp -rn $WORKSPACE_BUCKET/data/pgen_minimal_qc/plink_chr19_* .

# Get files for pheno/covar/step1
gsutil -m cp -rn $WORKSPACE_BUCKET/data/regenie_* .
gsutil -m cp -rn $WORKSPACE_BUCKET/data/regenie/* .

# Run regenie. I recommend the "--mcc" parameter for additional QC
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pfile plink_${curr_chr}_multi_split \
    --phenoFile regenie_pheno_zeroFID.txt \
    --covarFile regenie_covars_20pcs_withsex_zeroFID.txt \
    --bt --mcc \
    --firth --approx --pThresh 0.01 \
    --pred revised_20pcs_pred_sexcovar.list \
    --bsize 400 \
    --out aou_step2_rg_${curr_chr}_firthallvariants \
    --minMAC 20

# Backup results
gsutil -m cp -rn aou_step2_rg_${curr_chr}_firthallvariants $WORKSPACE_BUCKET/data/rg_results/
