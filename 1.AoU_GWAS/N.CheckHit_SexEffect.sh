# Operate on male or female sex to evaluate for role of sex on significance of hits
# Employs files produced by earlier steps
#    Plink file (QC is solely call rate and mac of 20): anc_hwe_test/${chr}

# Organize output
mkdir sex_eval

# Make file for keep function
awk '$4 == 0 {print $1 "\t" $2}' regenie_covar_20commonpcs.txt > female_ids.txt # Female is 0
awk '$4 == 1 {print $1 "\t" $2}' regenie_covar_20commonpcs.txt > male_ids.txt # Male is 1

for ((chr=1;chr<=22;chr++)); do \
  ./regenie_v3.2.8.gz_x86_64_Linux --step 2 --pgen anc_hwe_test/${chr} \
              --phenoFile regenie_pheno.txt --covarFile regenie_covar_20commonpcs.txt \
              --bt --firth-se --firth --approx --pThresh 0.01 --pred aou_step1_rg_array_anc_all_pred.list \
              --bsize 400 --minMAC 20 --phenoCol AD_any --keep female_ids.txt \
              --out sex_eval/${chr}_female ;\
  ./regenie_v3.2.8.gz_x86_64_Linux --step 2 --pgen anc_hwe_test/${chr} \
              --phenoFile regenie_pheno.txt --covarFile regenie_covar_20commonpcs.txt \
              --bt --firth-se --firth --approx --pThresh 0.01 --pred aou_step1_rg_array_anc_all_pred.list \
              --bsize 400 --minMAC 20 --phenoCol AD_any --keep male_ids.txt \
              --out sex_eval/${chr}_male ;\
done

# Organize data, produce GWAS for hits at out
head -n 1 sex_eval/${chr}_female_AD_any.regenie > sex_eval/aou_hits_female_stats.txt ;\
for file in sex_eval/*_female_AD_any.regenie; do \
    tail -n +2 "$file" >> sex_eval/aou_hits_female_stats.txt ;\
done

head -n 1 sex_eval/${chr}_male_AD_any.regenie > sex_eval/aou_hits_male_stats.txt ;\
for file in sex_eval/*_male_AD_any.regenie; do \
    tail -n +2 "$file" >> sex_eval/aou_hits_male_stats.txt ;\
done
