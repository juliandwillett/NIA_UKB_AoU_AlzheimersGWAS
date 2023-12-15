# Get files
gsutil -u $GOOGLE_PROJECT -m cp -r $WORKSPACE_BUCKET/data/array_data_for_regenie_step1/arrays_autosomes_post_qc_pruned* .
gsutil -m cp -rn $WORKSPACE_BUCKET/data/regenie_* .
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# Run regenie
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --qt --force-qt --mcc \
    --out aou_step1_rg_array \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20 \
    --phenoCol AD_any

# Run regenie with exclusion of related individuals
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --qt --force-qt --mcc \
    --out aou_step1_rg_array \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20 \
    --phenoCol AD_any
    --exclude relatedness_flagged_samples.tsv

# Backup results
gsutil -m cp -rn aou_step1_rg_array* $WORKSPACE_BUCKET/data/regenie_step1_mcc_qt/
