# I recommend doing this after exporting from AoU, given the limited file size you can successfully export after being approved for large egress.
head -n 1 aou_step2_rg_chr6_allvar_anc_all_AD_any_min20N_C.regenie > aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs.txt
for file in *.regenie; do \
    tail -n +2 "$file" >> aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs.txt ;\
done

##
# Do for ancestry stratified datasets
ancestries=(eur afr amr)
for anc in "${ancestries[@]}"; do \
    head -n 1 aou_step2_rg_chr7_allvar_anc_${anc}_AD_any_min20N_A.regenie > aou_AD_any_anc_${anc}_gwas.txt ;\
    for file in *${anc}*.regenie; do \
        tail -n +2 "$file" >> aou_AD_any_anc_${anc}_gwas.txt ;\
    done \
done
