# download files
mkdir gwas
gsutil -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/rg_results/* gwas/

outcomes=("AD" "Dementia" "Dementia_By_Proxy" "AD_any" "Dementia_any")

# merge files
for outcome in "${outcomes[@]}" ; do \
  echo "On outcome $outcome"
  head -n 1 gwas/aou_step2_rg_chr1_common_anc_all_${outcome}.regenie > aou_${outcome}_anc_all_gwas.txt ;\
  for ((i=1; i<=22; i++)); do \
    file="gwas/aou_step2_rg_chr${i}_common_anc_all_${outcome}.regenie" ;\
    tail -n +2 "$file" >> aou_${outcome}_anc_all_gwas.txt ;\
  done \
done

# convert to tsv, so it works with locus zoom
# Refer to Rmd file in FASRC

##################################

# run on ancestry stratified data in home computer directory
anc="eur"
for outcome in "${outcomes[@]}" ; do \
  echo "On outcome $outcome"  
  head -n 1 aou_step2_rg_chr1_common_anc_${anc}_${outcome}.regenie > aou_${outcome}_anc_${anc}_gwas.txt ;\
  for ((i=1; i<=22; i++)); do \
    file="aou_step2_rg_chr${i}_common_anc_${anc}_${outcome}.regenie" ;\
    tail -n +2 "$file" >> aou_${outcome}_anc_${anc}_gwas.txt ;\
  done \
done
