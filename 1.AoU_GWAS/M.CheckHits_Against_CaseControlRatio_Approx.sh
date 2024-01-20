# First ensure existing hits file is present
head -2 aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt

# Make folder to organize results
mkdir ccratio_approx_testing ; mkdir hwe_call_plink

# Make bed file for the hits
awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt > aou_hits.bed

# Randomly subsample controls to create file with case:control around 1:5 instead of around 1:20
awk 'NR==1 {print} $4 == 1 {print} $4 == 0 {if (rand() <= 0.2) print}' regenie_pheno.txt > regenie_pheno_cc_1to5.txt
awk '{count[$4]++} END {for (value in count) print value, count[value]}' regenie_pheno.txt

# Produce code to iterate and accomplish these tasks. First get QC'd data
for ((chr=1;chr<=22;chr++)) ; do \
  if (($chr>=17)) ; then
    ./plink2 --pfile plink_chr${chr}_multi_split \
        --geno 0.1 --hwe 1e-15 --set-all-var-ids @:#:\$r,\$a \
        --make-pgen --out hwe_call_plink/${chr} --new-id-max-allele-len 10000 \
        --extract bed1 aou_hits.bed ;\
  else \
    ./plink2 --pfile plink_chr${chr}_multi_split_merged \
        --geno 0.1 --hwe 1e-15 --set-all-var-ids @:#:\$r,\$a \
        --make-pgen --out hwe_call_plink/${chr} --new-id-max-allele-len 10000 \
        --extract bed1 aou_hits.bed ;\
  fi ;\
  if (($chr>=16)); then \ 
                 awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' hwe_call_plink/${chr}.psam > t ;\
                 mv t hwe_call_plink/${chr}.psam ; \
  else \
                awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $2 "\t" "NA"}' hwe_call_plink/${chr}.psam > t ;\
                 mv t hwe_call_plink/${chr}.psam ; \
  fi ;\
done

# Now run the regenie to do the tests
for ((chr=1;chr<=22;chr++)) ; do
  # approx, unchanged case:control
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen hwe_call_plink/${chr} \
    --phenoFile regenie_pheno.txt --covarFile regenie_covar_20commonpcs.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred aou_step1_rg_array_anc_all_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_approx_cc_1_20 \
    --minMAC 20 --phenoCol AD_any ;\

    # no approx, unchanged case:control
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen hwe_call_plink/${chr} \
    --phenoFile regenie_pheno.txt --covarFile regenie_covar_20commonpcs.txt \
    --bt --firth-se --firth --pThresh 0.01 \
    --pred aou_step1_rg_array_anc_all_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_approx_cc_1_20 \
    --minMAC 20 --phenoCol AD_any ;\

    # approx, case:control 1:5
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen hwe_call_plink/${chr} \
    --phenoFile regenie_pheno_cc_1to5.txt --covarFile regenie_covar_20commonpcs.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred aou_step1_rg_array_anc_all_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_approx_cc_1_20 \
    --minMAC 20 --phenoCol AD_any ;\

    # no approx, case:control 1:5
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen hwe_call_plink/${chr} \
    --phenoFile regenie_pheno_cc_1to5.txt --covarFile regenie_covar_20commonpcs.txt \
    --bt --firth-se --firth --pThresh 0.01 \
    --pred aou_step1_rg_array_anc_all_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_approx_cc_1_20 \
    --minMAC 20 --phenoCol AD_any ;\
done


