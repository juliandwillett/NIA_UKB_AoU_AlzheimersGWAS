# Less population stratification in single ancestry data, so less concern of bias.
# Given the number of hits in the multiancestry data, check that data against single-ancestry data, 
#   removing hits that did not pass HWE in single-ancestry data.
# So all we need are the HWE p values for single-ancestry cohorts

# Prepare bed file
mkdir hwe_testing ; mkdir hwe_testing/hardy_out/
ancestries=(eur afr amr sas eas mid)
awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_gwsig.txt > hwe_testing/aou_hits.bed

# Run the hardy calculations
for ((chr=1;chr<=22;chr++)); do \
  for anc in "${ancestries[@]}"; do \
    ./plink2 --pfile pgen_geno_1e-1_mac_20/chr${chr} --keep ancestries/${anc}_ids.txt \
      --missing --hardy midp --out hwe_testing/hardy_out/${chr}_${anc} --extract bed1 hwe_testing/aou_hits.bed ;\
  done \
done

# Then merge the output
for anc in "${ancestries[@]}"; do \
  head -n 1 hwe_testing/hardy_out/1_${anc}.hardy > hwe_testing/hardy_out/${anc}_hwe_stats.txt ;\
  for file in hwe_testing/hardy_out/*${anc}.hardy; do \
    tail -n +2 "$file" >> hwe_testing/hardy_out/${anc}_hwe_stats.txt ;\
  done \
done
