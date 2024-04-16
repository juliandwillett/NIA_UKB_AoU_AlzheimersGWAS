# FORMAT NIAGADS SO COMPATIBLE WITH OTHER DATASETS
nia_file="/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/selfALL25660_a_pc_5JPCs_noage_ss.Affection.Status.glm.logistic.full_nohwefilter.dn8.gz"
zcat $nia_file | awk 'NR==1 {$3 = "ID"; $5 = "ALLELE1"; $6 = "ALLELE0"; $7 = "A1FREQ"; $8 = "BETA"; $9 = "SE"; \
  $10 = "Pval"; print $0} NR>1 {$3 = $1 "-" $2 "-" $6 "-" $5; print $0}' >\
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/nia_all_for_meta_nohwe_qc.txt

# Run a script with the METAL parameters
sbatch 1A_RunMETAL.sh

# functions to clean up data for further processing and expedite analysis (intersections with single GWAS, for example)
mv METAANALYSIS1.TBL aou_ukb_nia_allvar_meta_analysis.TBL

# Add CHR, POS, and FAVOR_ID columns, to make sorting easier and enable later intersections (ie intersection with Bellenguez et al for comparison)
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0 "\tCHR\tPOS\tIDrev"} NR>1 {split($1, values, "-"); $16 = values[1]; \
  $17 = values[2]; $18 = $16 "-" $17 "-" toupper($3) "-" toupper($2); $20 = $16 "-" $17; print $0}' \
  aou_ukb_nia_allvar_meta_analysis.TBL > aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols.TBL

# Isolate GW significant hits to focus the analysis (and make R code work more efficiently)
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0} NR>1 && $10 <= 5e-8 && $16 != 23 {print $0}' \
  aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols.TBL > \
  aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols_gw_sig.TBL

#############
# R processing and make intersect reference file: in R
data = vroom("aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols_gw_sig.TBL") %>%
  select(-HetISq,-HetChiSq) %>% filter((MaxFreq - MinFreq) < 0.4,CHR < 23) %>% arrange(CHR,POS)
vroom_write(data,"aou_ukb_nia_allvar_meta_qc_sorted_gwsig.txt") 
vroom_write(data %>% select(MarkerName),"working/favor_hits_aou_vs_ukb.txt",col_names=F) # for favor 
vroom_write(data %>% mutate(CHRPOS = glue("{CHR}-{POS}")),"working/meta_hits_for_intersects.txt")

#############
# Intersect the meta significant hits with each GWAS to make getting p values more efficient. Start with AoU here
awk 'NR==1 {$16 = "CHRPOS"; print} NR>1 {$16 = $1 "-" $2; print}' \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals.txt > \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos.txt ;\
awk 'NR==FNR{arr[$18]; next} $16 in arr' working/meta_hits_for_intersects.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos.txt > \
  working/meta_hits_aou_intersect.txt

# Do the same for UKB
awk 'NR==1 {$15 = "CHRPOS"; print} NR>1 {$15 = $1 "-" $2; print}' \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id.regenie > \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos.regenie ;\
awk 'NR==FNR{arr[$18]; next} $15 in arr' working/meta_hits_for_intersects.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos.regenie > \
  working/meta_hits_ukb_intersect.txt

# Intersect with NIAGADS to check for mutual hits
awk 'NR==FNR{arr[$18]; next} $15 in arr' working/meta_hits_for_intersects.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos.txt > \
  working/meta_hits_niagads_intersect.txt

# Then make files for Manhattan to ensure efficiency
awk 'NR==1 {print $0} NR>1 && $16 < 23 && $7 - $6 < 0.4 && $10 < 1e-1 {print $0}' \
  /n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols.TBL > \
  /n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/working/meta_qc_for_manhattan.txt
