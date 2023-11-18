# I recommend doing this after exporting from AoU, given the limited file size you can successfully export after being approved for large egress.
head -n 1 aou_step2_rg_chr6_allvar_anc_all_AD_any_min20N_C.regenie > aou_AD_any_anc_all_gwas.txt

for file in *.regenie; do \
    tail -n +2 "$file" >> aou_AD_any_anc_all_gwas.txt ;\
done
