# Get files
gsutil -u $GOOGLE_PROJECT -m cp -r $WORKSPACE_BUCKET/data/array_data_for_regenie_step1/arrays_autosomes_post_qc_pruned* .
gsutil -m cp -rn $WORKSPACE_BUCKET/data/regenie_* .

# Run regenie
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned \
    --phenoFile regenie_pheno_revised.txt \
    --covarFile regenie_covar_revised.txt \
    --bt \
    --out aou_step1_rg_array \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20

# Backup results
gsutil -m cp -rn aou_step1_rg_array* $WORKSPACE_BUCKET/data/regenie/
