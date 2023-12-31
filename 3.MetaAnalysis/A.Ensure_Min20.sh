##### Double check that results correspond to minimum 20 individuals, per AoU requirements for making summary stats public
### AoU
awk 'NR==1 {print} NR>1 && $6 * $7 >= 20 && (1-$6) * $7 >= 20 {print}' \
/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse.txt > \
/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse_ensure20.txt

### UKB
awk 'NR==1 {print} NR>1 && $6 * $7 >= 20 && (1-$6) * $7 >= 20 {print}' \
/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids.txt > \
/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_ensure20.txt

### NIAGADS
awk 'NR==1 {print} NR>1 && $7 * $13 >= 20 && (1-$7) * $13 >= 20 {print}' \
/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS/NIAGADS_meta.tsv > \
/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_ensure20.txt
