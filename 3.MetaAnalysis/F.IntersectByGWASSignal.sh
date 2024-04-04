# INTERSECT BY CHRPOS, INTERSECT BY ID IN R WHERE IT IS EASIER TO BE PRECISE
awk 'NR==FNR{arr[$1]; next} $1 in arr' working/gwas_gwsig_for_intersect_ids.txt \
  /n/home09/jwillett/true_lab_storage/00_AoU/aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols.TBL > \
  working/gwas_hits_meta_intersect.txt ;\
awk 'NR==FNR{arr[$3]; next} $16 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos.txt > \
  working/gwas_hits_aou_intersect.txt ;\
awk 'NR==FNR{arr[$3]; next} $15 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos.regenie > \
  working/gwas_hits_ukb_intersect.txt ;\
awk 'NR==FNR{arr[$3]; next} $15 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos.txt > \
  working/gwas_hits_niagads_intersect.txt
