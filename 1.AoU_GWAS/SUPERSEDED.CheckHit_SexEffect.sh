DEPRECATED

# Operate on male or female sex to evaluate for role of sex on significance of hits
# Employs files produced by earlier steps
#    Plink file (QC is solely call rate and mac of 20): anc_hwe_test/${chr}

# Organize output
mkdir sex_eval ; mkdir sex_eval/out/

# Make file for keep function
awk '$4 == 0 {print $1 "\t" $2}' regenie_covar_20commonpcs.txt > female_ids.txt # Female is 0
awk '$4 == 1 {print $1 "\t" $2}' regenie_covar_20commonpcs.txt > male_ids.txt # Male is 1

for ((chr=1;chr<=22;chr++)); do \
  ./plink2 --pfile pgen_geno_1e-1_mac_20/chr${chr} \
      --make-pgen --out sex_eval/${chr} --extract bed1 all_study_hits.bed ;\
  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' sex_eval/${chr}.psam > t ;\
  mv t sex_eval/${chr}.psam ;\
  ./regenie_v3.2.8.gz_x86_64_Linux --step 2 --pgen sex_eval/${chr} \
              --phenoFile regenie_pheno.txt --covarFile regenie_covar_20pcs.txt \
              --bt --firth-se --firth --approx --pThresh 0.01 --pred revised_pred.list \
              --bsize 400 --minMAC 20 --phenoCol AD_any --keep female_ids.txt \
              --out sex_eval/out/${chr}_female ;\
  ./regenie_v3.2.8.gz_x86_64_Linux --step 2 --pgen sex_eval/${chr} \
              --phenoFile regenie_pheno.txt --covarFile regenie_covar_20pcs.txt \
              --bt --firth-se --firth --approx --pThresh 0.01 --pred revised_pred.list \
              --bsize 400 --minMAC 20 --phenoCol AD_any --keep male_ids.txt \
              --out sex_eval/out/${chr}_male ;\
done

# Organize data, produce GWAS for hits at out
sex=(female male) ;\
for s in "${sex[@]}"; do \
  head -n 1 sex_eval/out/1_${s}_AD_any.regenie > sex_eval/out/aou_hits_${s}_stats.txt ;\
  for file in sex_eval/out/*_${s}_AD_any.regenie; do \
      tail -n +2 "$file" >> sex_eval/out/aou_hits_${s}_stats.txt ;\
  done \
done
