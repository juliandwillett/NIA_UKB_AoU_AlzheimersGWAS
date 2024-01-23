# For working environment,  I recommend 32 CPUs with 120 GB RAM.

# Get necessary files
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip ;\
unzip plink2_linux*
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# Get files for current chromosome
bucket="gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67"
curr_chr="chr19"
gsutil -m cp -rn $bucket/data/pgen_minimal_qc/plink_${curr_chr}_* .

# Do QC, if it has not been run already. Did not do mind because AoU already flagged individuals and aimed to maximize sample size
# IMPORTANT: ENSURE THAT FID IID COLUMNS IN RELATED FLAGGED FOR REGENIE FILE MATCH PSAM FILE
./plink2 --pfile plink_${curr_chr}_multi_split_merged \
        --make-pgen --out plink_chr${chr}_allvar_anc_all \
        --mac 20

# Do extra bit of QC
./plink2 --pfile plink_out/plink_chr${chr}_allvar_anc_all \
        --make-pgen --out plink_out/plink_chr${chr}_allvar_anc_all_geno_1e-1 \
        --mac 20 --geno 0.1

# revise psam file given the empty column being dropped
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' plink_${curr_chr}_allvar_anc_all.psam > tmp
mv tmp plink_${curr_chr}_allvar_anc_all.psam

# Get files for pheno/covar/step1
gsutil -m cp -rn $bucket/data/regenie_* .
gsutil -m cp -rn $bucket/data/regenie/* .

# Rename workspace in pred file, to enable running on parallel workspaces
awk '{gsub("duplicateofalzheimersgwastake5", "duplicateofalzheimersgwastake5", $2)} 1' aou_step1_rg_array_anc_all_pred.list > revised_pred.list

# run regenie for non ancestry stratified in parallel
for ((chr=1;chr<=22;chr++)); do \
                ./regenie_v3.2.8.gz_x86_64_Linux \
                    --step 2 \
                    --pgen plink_out/plink_chr${chr}_allvar_anc_all_geno_1e-1 \
                    --phenoFile regenie_pheno.txt \
                    --covarFile regenie_covar_20commonpcs.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred aou_step1_rg_array_anc_all_pred.list \
                    --bsize 400 \
                    --out aou_rg_chr${chr}_firthallvar_geno_1e-1_commonpcs
                    --minMAC 20 \
                    --phenoCol AD_any ;\
                gsutil -m cp -r aou_step2_rg_${curr_chr}_firthallvariants_commonrarepcs* $WORKSPACE_BUCKET/data/rg_results_commonrarepcs/ ;\
done

# run regenie for ancestry-stratified cohorts
gsutil cp $WORKSPACE_BUCKET/data/*_ids.txt .
ancestries=(eur afr amr)
for ((chr=1;chr<=22;chr++)); do \
        curr_chr="chr${chr}" ;\
        for anc in "${ancestries[@]}"; do \
                awk '{gsub("duplicateofalzheimersgwastake5d1", "duplicateofalzheimersgwastake5", $2)} 1' aou_step1_rg_array_anc_${anc}_pred.list > revised_pred.list
                ./regenie_v3.2.8.gz_x86_64_Linux \
                    --step 2 \
                    --pgen plink_${curr_chr}_allvar_anc_all \
                    --phenoFile regenie_pheno.txt \
                    --covarFile regenie_covar.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred revised_pred.list \
                    --bsize 400 \
                    --out aou_step2_rg_${curr_chr}_firthallvariants_${anc} \
                    --minMAC 20 \
                    --phenoCol AD_any --keep ${anc}_ids.txt ;\
                gsutil -m cp -r aou_step2_rg_${curr_chr}_firthallvariants_${anc}* $WORKSPACE_BUCKET/data/rg_results_rel_anc_strat/ ;\
        done  
done

#####
# Code for parallelizing
