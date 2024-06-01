# Random case control assignment to check for noise, 5% case ratio

{R}
tmp = vroom('regenie_input/regenie_pheno.txt')
tmp %<>% mutate(AD_any = ifelse(runif(n()) < 0.05, 1, 0))
head(tmp)
table(tmp$AD_any)
vroom_write(tmp,'regenie_input/regenie_pheno_randomcases_5percent.txt')
{/R}

# Step 1
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/arrays_autosomes_post_qc_pruned_common \
    --phenoFile regenie_input/regenie_pheno_randomcases_5percent.txt \
    --covarFile regenie_input/regenie_covar_commonpcs.txt \
    --bt \
    --out step1_randomcase/aou_step1_randomcase \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_40 \
    --phenoCol AD_any
