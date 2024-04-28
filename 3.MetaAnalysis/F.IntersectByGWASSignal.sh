# INTERSECT BY CHRPOS FOR SPEED FIRST, THEN INTERSECT BY ID IN R WHERE IT IS EASIER TO BE PRECISE
# NIAGADS INTERSECT
awk 'NR==FNR{arr[$3]; next} $15 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_chrpos.txt > \
  working/niagads_intersects_chrpos.txt ;\

# NIMH INTERSECT:
awk 'NR==FNR{arr[$3]; next} $15 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIMH/minorallele_merged_res_NIMH+NIA2247_Aff.offset.a.full.dn8_formeta_chrpos.txt > \
  working/nimh_intersects_chrpos.txt ;\

# NIAGADS NIMH META INTERSECT
awk 'NR==FNR{arr[$3]; next} $19 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/00_AoU/final_meta_results/meta_niagads_nimh_chrposrefalt_chrpos.TBL > \
  working/niagads_nimh_meta_intersects_chrpos.txt ;\

# UKB intersect
awk 'NR==FNR{arr[$3]; next} $15 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos.regenie > \
  working/ukb_intersects_chrpos.txt ;\

# AOU INTERSECT:
awk 'NR==FNR{arr[$3]; next} $15 in arr' working/gwas_gwsig_for_intersect.txt \
  ../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_chrpos.txt > \
  working/aou_intersects_chrpos.txt ;\

# UKB AOU META INTERSECT
awk 'NR==FNR{arr[$3]; next} $19 in arr' working/gwas_gwsig_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/00_AoU/final_meta_results/meta_ukb_aou_chrposrefalt_chrpos.TBL > \
  working/ukb_aou_meta_intersects_chrpos.txt ;\

# NIAGADS NIMH UKB AOU META INTERSECT
awk 'NR==FNR{arr[$3]; next} $19 in arr' working/gwas_gwsig_for_intersect.txt \
  final_meta_results/meta_niagads_nimh_ukb_aou_chrposrefalt_chrpos.TBL > \
  working/ukb_aou_meta_intersects_chrpos.txt
