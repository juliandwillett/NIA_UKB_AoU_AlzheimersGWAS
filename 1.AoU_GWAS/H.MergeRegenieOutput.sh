# I recommend doing this after exporting from AoU, given the limited file size you can successfully export after being approved for large egress.
head -n 1 aou_step2_rg_chr16_AD_any.regenie > aou_AD_any_anc_all_gwas.txt
for ((i=1; i<=22; i++)); do \
    file="aou_step2_rg_chr${i}_AD_any.regenie" ;\
    tail -n +2 "$file" >> aou_AD_any_anc_all_gwas.txt ;\
done
gsutil cp -rn aou_AD_any_anc_all_gwas.txt $WORKSPACE_BUCKET/data/rg_results/
