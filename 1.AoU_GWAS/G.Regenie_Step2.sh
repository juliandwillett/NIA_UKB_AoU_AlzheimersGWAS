# For working environment,  I recommend 32 CPUs with 120 GB RAM.

# Get necessary files
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip ;\
unzip plink2_linux*
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# Make requisite folders
mkdir pgen_files ; mkdir pgen_geno_1e-1_mac_20

# Get files for current chromosome
bucket="gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67"
gsutil -m cp -rn $bucket/data/pgen_minimal_qc/plink_* pgen_files/

# Do QC. Did not do mind because AoU already flagged individuals and aimed to maximize sample size
for ((chr=4;chr<=22;chr++)); do \
  if (($chr>=17)) ; then \
      ./plink2 --pfile pgen_files/plink_chr${chr}_multi_split \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 --memory 200000 \
          --make-pgen --out pgen_geno_1e-1_mac_20/chr${chr} --new-id-max-allele-len 10000 ;\
    else \
      ./plink2 --pfile pgen_files/plink_chr${chr}_multi_split_merged \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 --memory 200000 \
          --make-pgen --out pgen_geno_1e-1_mac_20/chr${chr} --new-id-max-allele-len 10000 ;\
  fi ;\
  if (($chr>=16)); then \ 
                   awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' pgen_geno_1e-1_mac_20/chr${chr}.psam > t ;\
                   mv t pgen_geno_1e-1_mac_20/chr${chr}.psam ; \
    else \
                  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $2 "\t" "NA"}' pgen_geno_1e-1_mac_20/chr${chr}.psam > t ;\
                   mv t pgen_geno_1e-1_mac_20/chr${chr}.psam ; \
  fi ;\
done

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
