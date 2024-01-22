# Less population stratification in single ancestry data, so less concern of bias.
# Given the number of hits in the multiancestry data, check that data against single-ancestry data, 
#   removing hits that did not pass HWE in single-ancestry data.
# So all we need are the HWE p values for single-ancestry cohorts

# Prepare bed file
ancestries=(eur afr amr)
awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt > aou_hits.bed
mkdir anc_hwe_test ; mkdir anc_hwe_test/hardy_out/

for ((chr=1;chr<=22;chr++)); do \
  for anc in "${ancestries[@]}"; do \
     # Loop through chromosomes to do QC (just call rate)
     if (($chr>=17)) ; then \
      ./plink2 --pfile plink_chr${chr}_multi_split \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 \
          --make-pgen --out anc_hwe_test/${chr}_${anc} --new-id-max-allele-len 10000 \
          --extract bed1 aou_hits.bed --keep ${anc}_ids.txt ;\
    else \
      ./plink2 --pfile plink_chr${chr}_multi_split_merged \
          --geno 0.1 --set-all-var-ids @-#-\$r-\$a --mac 20 \
          --make-pgen --out anc_hwe_test/${chr}_${anc} --new-id-max-allele-len 10000 \
          --extract bed1 aou_hits.bed --keep ${anc}_ids.txt ;\
    fi ;\
    if (($chr>=16)); then \ 
                   awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' anc_hwe_test/${chr}_${anc}.psam > t ;\
                   mv t anc_hwe_test/${chr}_${anc}.psam ; \
    else \
                  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $2 "\t" "NA"}' anc_hwe_test/${chr}_${anc}.psam > t ;\
                   mv t anc_hwe_test/${chr}_${anc}.psam ; \
    fi ;\

    # Next run the hardy testing
    ./plink2 --pfile anc_hwe_test/${chr}_${anc} \
      --missing --hardy --out anc_hwe_test/hardy_out/${chr}_${anc} ;\
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
