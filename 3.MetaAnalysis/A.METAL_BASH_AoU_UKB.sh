# first ensure the ID columns are the same format in GWAS files, and create a P value column
awk -v OFS='\t' 'NR==1 {print $0 "\tPval"} NR>1 {$15 = 10^(-1 * $12); print }' /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs.txt > /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals.txt
awk -v OFS='\t' 'NR>1 {$3 = $1 ":" $2 ":" $4 "," $5}1' aou_AD_any_anc_all_gwas_pvals.txt > aou_AD_any_anc_all_gwas_pvals_ids.txt # split multiallelics have "." as ID, so fix this

# Add p val, IDs column to UKB
awk -v OFS='\t' 'NR==1 { $14 = "Pval"; print} NR>1 {$3 = $1 ":" $2 ":" $4 "," $5; $14 = 10^(-1 * $12); print }' \
/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_ensure20.txt > \
/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_pval_ensure20.txt

# Run a script with the METAL parameters
salloc -p test --mem 80000 -t 0-02:00 -n 4
metal METAL_script.txt

# functions to clean up data for further processing and expedite analysis (intersections with single GWAS, for example)
mv METAANALYSIS1.TBL aou_ukb_allvar_meta_analysis.TBL

# Add CHR, POS, and FAVOR_ID columns, to make sorting easier and enable later intersections (ie intersection with Bellenguez et al for comparison)
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0 "\tCHR\tPOS\tID\tIDrev"} NR>1 {split($1, values, ":"); $16 = values[1]; $17 = values[2]; $18 = $16 "-" $17 "-" toupper($2) "-" toupper($3); $19 = $16 "-" $17 "-" toupper($3) "-" toupper($2); $20 = $16 "-" $17; print $0}' \
  aou_ukb_allvar_meta_analysis.TBL > aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL

# Isolate GW significant hits to focus the analysis (and make R code work more efficiently)
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0} NR>1 && $10 <= 5e-8 {print $0}' aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL > \
  aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols_gw_sig.TBL

# to make intersect reference file: in R
data = vroom("aou_ukb_allvar_meta_qc_sorted_gwsig.txt")
data %<>% mutate(CHRPOS = glue("{CHR}-{POS}"))
vroom_write(data,"meta_hits_for_intersects_aou_vs_ukb.txt")

# Intersect the meta significant hits with each GWAS to make getting p values more efficient
awk 'NR==FNR{arr[$19]; next} $16 in arr' meta_hits_for_intersects_aou_vs_ukb.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse_ensure20.txt > \
  meta_hits_aou_intersect_aou_vs_ukb.txt
awk 'NR==1 {$14 = "CHROMPOS"; print} NR>1 {$14 = $1 "-" $2; print}' /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_ensure20.txt > \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_ensure20_chrpos.txt
awk 'NR==FNR{arr[$19]; next} $14 in arr' meta_hits_for_intersects_aou_vs_ukb.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_ensure20_chrpos.txt > \
  meta_hits_ukb_intersect_aou_vs_ukb.txt

# Intersect with NIAGADS to check for mutual hits
awk 'NR==FNR{arr[$19]; next} $15 in arr' meta_hits_for_intersects_aou_vs_ukb.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos.txt > \
  meta_hits_niagads_intersect_aou_vs_ukb.txt

# Also intersect with Bellenguez sumstats to speed up checking for mutual hits, re concern for FP
awk 'NR==FNR{arr[$19]; next} $18 in arr' meta_hits_for_intersects_aou_vs_ukb.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/Other_GWAS/bellengeuz_meta_ids_chrpos.tsv > \
  meta_hits_bell_intersect_aou_vs_ukb.txt

# Then make files for Manhattan to ensure efficiency
awk 'NR==1 {print $0} NR>1 && $16 < 23 && $7 - $6 < 0.4 && $10 < 1e-3 {print $0}' \
/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL > \
/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/meta_qc_for_manhattan.txt
