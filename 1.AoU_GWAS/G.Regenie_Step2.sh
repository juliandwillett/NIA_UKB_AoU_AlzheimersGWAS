# For plink work,  I recommend 32 CPUs with 120 GB RAM.
# For regenie step 1, 16 cores 104 GB.
# For regenie work, I recommend 64 CPUs with 204 GB of RAM. 
#  32 cores with matching ram is 50% less cost with 50% less speed, so 64 makes more sense.
#  96 cores with 86 ram is about the same speed as 64 with 204.
# 400 GB of RAM did not make any difference vs 200 GB

# Get necessary files
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip ;\
unzip plink2_linux*
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# Make requisite folders
mkdir pgen_files ; mkdir pgen_geno_1e-1_mac_20 ; mkdir rg_multi_geno_1e-1_mac_20
mkdir rg_multi_geno_1e-1_mac_20_for_export

# Get files for current chromosome
bucket="gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67"
gsutil -m cp -rn $bucket/data/pgen_minimal_qc/plink_* pgen_files/

# Do QC. Did not do mind because AoU already flagged individuals and aimed to maximize sample size
for ((chr=5;chr<=22;chr++)); do \
  if (($chr>=17)) ; then \
      ./plink2 --pfile pgen_files/plink_chr${chr}_multi_split \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 --memory 300000 \
          --make-pgen --out pgen_geno_1e-1_mac_20/chr${chr} --new-id-max-allele-len 10000 ;\
    else \
      ./plink2 --pfile pgen_files/plink_chr${chr}_multi_split_merged \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 --memory 300000 \
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

# run regenie for non ancestry stratified in parallel. While not QC'd on HWE for step 2, HWE QC was done for step 1 to get population structure 
for ((chr=3;chr<=22;chr++)); do \
                ./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/chr${chr}_geno_mac \
                    --phenoFile regenie_input/regenie_pheno.txt \
                    --covarFile regenie_input/regenie_covar_commonpcs_with_anc.txt \
                    --bt --firth-se --catCovarList ancestry_pred \
                    --firth --approx --pThresh 0.01 \
                    --pred aou_step1_rg_array_commonpcs_anccovarcat/_pred.list \
                    --bsize 400 \
                    --out aou_step2_rg_array_commonpcs_anccovarcat/chr${chr} \
                    --minMAC 20 \
                    --phenoCol AD_any ;\
        gsutil -m cp -r rg_multi_geno_1e-1_mac_20/chr${chr}* $bucket/data/rg_results_20commonpcs_geno_1e-1_mac_20_no_hwe/ ;\
done

### run regenie for ancestry-stratified cohorts
# produce files with FID/IID for each ancestry (determined by AoU)
mkdir ancestries
gsutil -u $GOOGLE_PROJECT -m cp -rn gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv ancestries/
ancestries=(eur afr amr sas eas mid) ;\
for anc in "${ancestries[@]}"; do \
  awk -v a="$anc" '$2 == a {print "0\t" $1 "\t" $2}' ancestries/ancestry_preds.tsv > ancestries/${anc}_ids.txt ;\
done
gsutil -m cp -rn ancestries/*ids.txt $bucket/data/ancestry_ids/

# get backed up data and run
gsutil cp $WORKSPACE_BUCKET/data/*_ids.txt .
ancestries=(amr afr eur)
for ((chr=2;chr<=22;chr++)); do \
        curr_chr="chr${chr}" ;\
        for anc in "${ancestries[@]}"; do \
                # orig string in first line, replacement in second
                #awk '{gsub("duplicateofalzheimersgwastake5d1", "duplicateofduplicateofalzheimersgwastake6", $2)} 1' aou_step1_rg_array_anc_${anc}_pred.list > revised_pred.list
                ./regenie_v3.2.8.gz_x86_64_Linux \
                    --step 2 \
                    --pgen pgen_geno_1e-1_mac_20/chr${chr} \
                    --phenoFile regenie_pheno.txt \
                    --covarFile regenie_covar_${anc}.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_step1_singleanc/aou_step1_rg_common_${anc}_pred.list \
                    --bsize 400 \
                    --out rg_step2_singleanc_anc_specific_pcs/aou_step2_rg_${curr_chr}_firth_${anc}_ancspecific_pcs \
                    --minMAC 20 \
                    --phenoCol AD,AD_any --keep ancestries/${anc}_ids.txt ;\
        done  
done
gsutil -m cp -r rg_step2_singleanc_anc_specific_pcs/* $WORKSPACE_BUCKET/data/rg_results_anc_strat_ancspecific_pcs/ ;\

# Check single variants for sex specific effects
./regenie_v3.4.1.gz_x86_64_Centos7_mkl \
                    --step 2 \
                    --pgen pgen_qc/18_48238744 \
                    --phenoFile regenie_input/regenie_pheno.txt \
                    --covarFile regenie_input/regenie_covar_commonpcs_female.txt \
                    --bt --firth-se \
                    --firth --approx --pThresh 0.01 \
                    --pred rg_sexspecific/rg_step1_female_pred.list \
                    --bsize 400 \
                    --out rg_sexspecific/rg_step2_chr18_48238744_female \
                    --minMAC 20 \
                    --phenoCol AD_any 
