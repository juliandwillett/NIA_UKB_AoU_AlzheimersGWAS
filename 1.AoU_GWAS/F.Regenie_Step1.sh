# Get files
gsutil -u $GOOGLE_PROJECT -m cp -r $WORKSPACE_BUCKET/data/array_data_for_regenie_step1/arrays_autosomes_post_qc_pruned_common* .
gsutil -m cp -rn $WORKSPACE_BUCKET/data/regenie_* .
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# Add header to psam file if needed
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' arrays_autosomes_post_qc_pruned_common.psam > tmp
mv tmp arrays_autosomes_post_qc_pruned_common.psam

# Run regenie
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned_common \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt \
    --out aou_step1_rg_array_common_rare_pcs \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any

# Run regenie with exclusion of related individuals (results less relevant when using regenie, that accounts for relatedness to an extent)
awk 'NR==1 {print "#FID\tIID"} NR>1 {print "0\t" $1}' relatedness_flagged_samples.tsv > related_flagged_for_regenie.txt
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt \
    --out aou_step1_rg_array_norelated \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20 \
    --phenoCol AD_any \
    --remove related_flagged_for_regenie.txt

# Rerun regenie step 1 by ancestry, using appropriate covar file (PCs for single ancestry)
ancestries=(eur amr afr) ;\
for anc in "${ancestries[@]}"; do \
    ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen piezo2_work/array_data/arrays_autosomes_post_qc_pruned_common \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar_${anc}.txt \
    --bt \
    --out rg_step1_singleanc/aou_step1_rg_common_${anc} \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20 \
    --phenoCol AD,AD_any ;\
done
gsutil -m cp -rn rg_step1_singleanc/aou_step1_rg_common* $WORKSPACE_BUCKET/data/regenie_step1_singleanc_anc_pcs/

# Backup results
gsutil -m cp -rn aou_step1_rg_array_norelated_* $WORKSPACE_BUCKET/data/regenie_step1_norelated/

