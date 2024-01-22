# Less population stratification in single ancestry data, so less concern of bias.
# Given the number of hits in the multiancestry data, check that data against single-ancestry data, 
#   removing hits that did not pass HWE in single-ancestry data.
# So all we need are the HWE p values for single-ancestry cohorts

# Prepare bed file
ancestries=(eur afr amr)
awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt > aou_hits.bed
mkdir anc_hwe_test ; mkdir anc_hwe_test/hardy_out/

for ((chr=1;chr<=22;chr++)); do \
  if (($chr>=17)) ; then \
      ./plink2 --pfile plink_chr${chr}_multi_split \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 \
          --make-pgen --out anc_hwe_test/${chr} --new-id-max-allele-len 10000 \
          --extract bed1 aou_hits.bed ;\
    else \
      ./plink2 --pfile plink_chr${chr}_multi_split_merged \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 \
          --make-pgen --out anc_hwe_test/${chr} --new-id-max-allele-len 10000 \
          --extract bed1 aou_hits.bed ;\
  fi ;\
  if (($chr>=16)); then \ 
                   awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' anc_hwe_test/${chr}.psam > t ;\
                   mv t anc_hwe_test/${chr}.psam ; \
    else \
                  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $2 "\t" "NA"}' anc_hwe_test/${chr}.psam > t ;\
                   mv t anc_hwe_test/${chr}.psam ; \
  fi ;\

  # now do the HWE by ancestry
  for anc in "${ancestries[@]}"; do \
    ./plink2 --pfile anc_hwe_test/${chr} --keep ${anc}_ids.txt \
      --missing --hardy midp --out anc_hwe_test/hardy_out/${chr}_${anc} ;\
  done \
done

# Then merge the output
# eur
head -n 1 anc_hwe_test/hardy_out/1_eur.hardy > anc_hwe_test/hardy_out/eur_hwe_stats.txt ;\
for file in anc_hwe_test/hardy_out/*eur.hardy; do \
    tail -n +2 "$file" >> anc_hwe_test/hardy_out/eur_hwe_stats.txt ;\
done
# afr
head -n 1 anc_hwe_test/hardy_out/1_afr.hardy > anc_hwe_test/hardy_out/afr_hwe_stats.txt ;\
for file in anc_hwe_test/hardy_out/*afr.hardy; do \
    tail -n +2 "$file" >> anc_hwe_test/hardy_out/afr_hwe_stats.txt ;\
done
# amr
head -n 1 anc_hwe_test/hardy_out/1_amr.hardy > anc_hwe_test/hardy_out/amr_hwe_stats.txt ;\
for file in anc_hwe_test/hardy_out/*amr.hardy; do \
    tail -n +2 "$file" >> anc_hwe_test/hardy_out/amr_hwe_stats.txt ;\
done
