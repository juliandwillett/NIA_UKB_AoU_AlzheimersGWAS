DEPRECATED

################
# GET FOR ANC STRATIFIED RESULTS
ancestries=(eur amr afr) ;\

# Only autosomes for step 1, do chr 23 on step 2
for anc in "${ancestries[@]}"; do \
  ./plink2 --pfile array_data/arrays_allchr --maf 0.01 --mac 100 --geno 0.1 --keep regenie_input/${anc}_ids.txt \
    --make-pgen --chr 1-22 --out array_data/arrays_${anc}_autosomes_post_qc --indep-pairwise 100kb 1 0.1 ;\
  ./plink2 --pfile array_data/arrays_${anc}_autosomes_post_qc --exclude array_data/arrays_${anc}_autosomes_post_qc.prune.out \
    --make-pgen --out array_data/arrays_${anc}_autosomes_post_qc_pruned ;\
  ./plink2 --pfile array_data/arrays_${anc}_autosomes_post_qc_pruned --pca 20 approx \
    --out array_data/arrays_${anc}_autosomes_post_qc_pruned_pcs ;\
done

# R WORK
covar = vroom("regenie_input/regenie_covar.txt",show_col_types = F) %>% select(FID,IID,Age,Sex)
for (anc in c("eur", "amr", "afr")) {
    pcs = vroom(glue("array_data/arrays_{anc}_autosomes_post_qc_pruned_pcs.eigenvec")) %>% rename(IID = `#IID`)
    anc_covar = merge(covar,pcs,by='IID') %>% relocate(FID,.before='IID')
    vroom_write(anc_covar,glue("regenie_input/regenie_covar_{anc}.txt"))
}

# REGENIE STEP 1
ancestries=(eur amr afr) ;\
for anc in "${ancestries[@]}"; do \
  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' array_data/arrays_${anc}_autosomes_post_qc_pruned.psam > tmp ;\
  mv tmp array_data/arrays_${anc}_autosomes_post_qc_pruned.psam ;\
done

for anc in "${ancestries[@]}"; do \
    ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
    --step 1 \
    --pgen array_data/arrays_${anc}_autosomes_post_qc_pruned \
    --phenoFile regenie_input/regenie_pheno.txt \
    --covarFile regenie_input/regenie_covar_${anc}.txt \
    --bt \
    --out rg_step1_singleanc/aou_step1_rg_common_${anc} \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_rg_20 \
    --phenoCol AD_any ;\
done ;\
gsutil -m cp -rn rg_step1_singleanc/aou_step1_rg_common* $WORKSPACE_BUCKET/data/regenie_step1_singleanc_anc_pcs/

# REGENIE STEP 2
ancestries=(amr afr eur)
for ((chr=23;chr<=23;chr++)); do \
        curr_chr="chr${chr}" ;\
        for anc in "${ancestries[@]}"; do \
                # orig string in first line, replacement in second
                #awk '{gsub("duplicateofalzheimersgwastake5d1", "duplicateofduplicateofalzheimersgwastake6", $2)} 1' aou_step1_rg_array_anc_${anc}_pred.list > revised_pred.list
                ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/chr${chr}_geno_mac \
                    --phenoFile regenie_input/regenie_pheno.txt \
                    --covarFile regenie_input/regenie_covar_${anc}.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_${anc}_pred.list \
                    --bsize 400 \
                    --out rg_step2_singleanc/${curr_chr}_${anc}_ancspecific_pcs \
                    --minMAC 20 \
                    --phenoCol AD_any ;\
        done  ;\
done
gsutil -m cp -r rg_step2_singleanc/* $WORKSPACE_BUCKET/data/regenie_step2_singleanc_anc_pcs/ ;\
