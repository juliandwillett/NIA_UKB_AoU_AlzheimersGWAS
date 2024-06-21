DEPRECATED

# First ensure existing hits file is present
head -2 aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt

# Make folder to organize results
mkdir ccratio_approx_testing ; mkdir hwe_call_plink

# Make bed file for the hits
#awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt > aou_hits.bed

# Randomly subsample controls to create file with case:control around 1:5 instead of around 1:20
awk 'NR==1 {print} $4 == 1 {print} $4 == 0 {if (rand() <= 0.2) print}' regenie_pheno.txt > regenie_pheno_cc_1to5.txt ;\
awk '{count[$4]++} END {for (value in count) print value, count[value]}' regenie_pheno.txt

# Produce code to iterate and accomplish these tasks. Use QC'd data produced by prior steps:
# sex_eval/${chr}

# Now run regenie to do the tests
for ((chr=1;chr<=22;chr++)) ; do
  # approx, unchanged case:control
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen sex_eval/${chr} \
    --phenoFile regenie_pheno.txt --covarFile regenie_covar_20pcs.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred revised_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_approx_cc_1_20 \
    --minMAC 20 --phenoCol AD_any ;\

    # no approx, unchanged case:control
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen sex_eval/${chr} \
    --phenoFile regenie_pheno.txt --covarFile regenie_covar_20pcs.txt \
    --bt --firth-se --firth --pThresh 0.01 \
    --pred revised_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_notapprox_cc_1_20 \
    --minMAC 20 --phenoCol AD_any ;\

    # approx, case:control 1:5
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen sex_eval/${chr} \
    --phenoFile regenie_pheno_cc_1to5.txt --covarFile regenie_covar_20pcs.txt \
    --bt --firth-se --firth --approx --pThresh 0.01 \
    --pred revised_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_approx_cc_1_5 \
    --minMAC 20 --phenoCol AD_any ;\

    # no approx, case:control 1:5
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 --pgen sex_eval/${chr} \
    --phenoFile regenie_pheno_cc_1to5.txt --covarFile regenie_covar_20pcs.txt \
    --bt --firth-se --firth --pThresh 0.01 \
    --pred revised_pred.list --bsize 400 \
    --out ccratio_approx_testing/chr${chr}_notapprox_cc_1_5 \
    --minMAC 20 --phenoCol AD_any ;\
done

names=(approx_cc_1_20 notapprox_cc_1_20 approx_cc_1_5 notapprox_cc_1_5) ;\
for name in "${names[@]}"; do \
  echo $name ;\
  head -n 1 ccratio_approx_testing/chr1_${name}_AD_any.regenie > ccratio_approx_testing/aou_hits_${name}.txt ;\
  for file in ccratio_approx_testing/*${name}_AD_any.regenie; do \
      tail -n +2 "$file" >> ccratio_approx_testing/aou_hits_${name}.txt ;\
  done \
done

###########
# R code to organize the results
aou_hits = vroom("aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_gwsig.txt",show_col_types=F) %>% 
    arrange(CHROM,GENPOS)
approx_cc_df = data.frame(ID=aou_hits$ID,P_Approx_120=NA,P_NotApprox_120=NA,P_Approx_15=NA,P_NotApprox_15=NA) 

options = c("approx_cc_1_20","notapprox_cc_1_20","approx_cc_1_5","notapprox_cc_1_5")
approx_data = lapply(X = options,
                     FUN = function(x) {
                         vroom(glue("ccratio_approx_testing/aou_hits_{x}.txt"),show_col_types=F) %>%
                             arrange(CHROM,GENPOS) %>% distinct(ID,.keep_all=T)
                     })

for (row in 1:nrow(approx_cc_df)) {
    curr_id = approx_cc_df$ID[[row]]    
    for (op in 1:4) {
        curr_opt = which(approx_data[[op]]$ID == curr_id)
        if (length(curr_opt) > 0) 
            approx_cc_df[[row,1+op]] = 10^(-approx_data[[op]]$LOG10P[[curr_opt]])
    }
}
head(approx_cc_df)
vroom_write(approx_cc_df,"aou_hits_approx_vs_cc_1_20_v_1_5_effect.txt")
