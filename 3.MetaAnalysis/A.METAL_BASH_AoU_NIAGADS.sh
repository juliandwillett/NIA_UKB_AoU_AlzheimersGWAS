# Run a script with the METAL parameters
salloc -p test --mem 80000 -t 0-02:00 -n 4
metal METAL_script_AOU_NIAGADS.txt

# functions to clean up data for further processing and expedite analysis (intersections with single GWAS, for example)
mv METAANALYSIS1.TBL aou_niagads_allvar_meta_analysis.TBL

# Add CHR, POS, and FAVOR_ID columns, to make sorting easier and enable later intersections (ie intersection with Bellenguez et al for comparison)
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0 "\tCHR\tPOS\tID\tIDrev"} NR>1 {split($1, values, ":"); $16 = values[1]; $17 = values[2]; $18 = $16 "-" $17 "-" toupper($2) "-" toupper($3); $19 = $16 "-" $17 "-" toupper($3) "-" toupper($2); $20 = $16 "-" $17; print $0}' \
  aou_niagads_allvar_meta_analysis.TBL > aou_niagads_allvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL

# Isolate GW significant hits to focus the analysis (and make R code work more efficiently)
awk 'BEGIN{FS=" "; OFS="\t"} NR==1 {print $0} NR>1 && $10 <= 5e-8 {print $0}' aou_niagads_allvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL > \
  aou_niagads_allvar_meta_analysis_IDcolon_chrposrefalt_cols_gw_sig.TBL

# to make intersect reference file: in R
data = vroom("aou_niagads_allvar_meta_qc_sorted_gwsig.txt") # make in R using first block of code
data %<>% mutate(CHRPOS = glue("{CHR}-{POS}"))
vroom_write(data,"meta_hits_for_intersects_vsniagads.txt")

# Intersect the meta significant hits with each GWAS to make getting p values more efficient
awk 'NR==FNR{arr[$19]; next} $16 in arr' meta_hits_for_intersects_vsniagads.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse_ensure20.txt > \
  meta_hits_aou_intersect_vs_niagads.txt
awk 'NR==1 {$15 = "CHRPOS"; print} NR>1 {$15 = $1 "-" $2; print}' /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_ensure20.txt > \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos_ensure20.txt
awk 'NR==FNR{arr[$19]; next} $15 in arr' meta_hits_for_intersects_vsniagads.txt \
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos_ensure20.txt > \
  meta_hits_niagads_intersect_vs_niagads.txt

# Then make files for Manhattan to ensure efficiency
awk 'NR==1 {print $0} NR>1 && $16 < 23 && $7 - $6 < 0.4 && $10 < 1e-3 {print $0}' \
/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL > \
/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/meta_qc_for_manhattan.txt
