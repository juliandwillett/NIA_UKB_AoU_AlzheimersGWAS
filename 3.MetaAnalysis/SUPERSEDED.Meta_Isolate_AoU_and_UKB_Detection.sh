{R} 
meta = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols_p_0_01.TBL",
  show_col_types = F) %>% mutate(CHRPOS = glue("{CHR}-{POS}")) %>% select(CHRPOS)
vroom_write(meta_chrpos,"meta_chr_pos_for_intersect.txt")
{/R}

# now intersect these chrpos against AoU and UKB
awk 'NR==FNR{arr[$1]; next} $16 in arr' meta_chr_pos_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos.txt > \
  meta_aou_intersect_in_aou_and_ukb.txt

# Do the same for UKB
awk 'NR==FNR{arr[$1]; next} $15 in arr' meta_chr_pos_for_intersect.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos.regenie > \
  meta_ukb_intersect_in_aou_and_ukb.txt

# Then plot in R
