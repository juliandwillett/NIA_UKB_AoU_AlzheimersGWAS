# download files
mkdir gwas
gsutil -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/rg_results/* gwas/

outcomes=("AD" "Dementia" "Dementia_By_Proxy" "AD_any" "Dementia_any")

# merge files
for outcome in "${outcomes[@]}" ; do \
  echo "On outcome $outcome"
  head -n 1 gwas/aou_step2_rg_chr1_common_anc_all_AD_any.regenie > aou_AD_any_anc_all_gwas.txt ;\
  for ((i=1; i<=22; i++)); do \
    file="gwas/aou_step2_rg_chr${i}_common_anc_all_AD_any.regenie" ;\
    tail -n +2 "$file" >> aou_AD_any_anc_all_gwas.txt ;\
  done \
done


