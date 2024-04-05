# Less population stratification in single ancestry data, so less concern of bias.
# Given the number of hits in the multiancestry data, check that data against single-ancestry data, 
#   removing hits that did not pass HWE in single-ancestry data.
# So all we need are the HWE p values for single-ancestry cohorts

###########
# FASRC: make files with p <= 1e-5 for Manhattans and QC
awk 'NR==1 || $10 <= 1e-5 {print}' /n/home09/jwillett/true_lab_storage/00_AoU/aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols.TBL >\
  /n/home09/jwillett/true_lab_storage/00_AoU/aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols_p_1e-5.TBL ;\
awk 'NR==1 || $14 <= 1e-5 {print}' /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos.regenie >\
  /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_p_1e-5.regenie ;\
awk 'NR==1 || $15 <= 1e-5 {print}' /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos.txt >\
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos_p_1e-5.txt ;\
awk 'NR==1 || $10 <= 1e-5 {print}' /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos.txt >\
  /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos_p_1e-5.txt

###########
#{R} FASRC
make_files_for_hwe = function() {
  # get CHRPOS of key variants for HWE testing. Not using IDS here as it is slower
  # The output also includes IDs, so I can do matching there.
  meta_hits = vroom("aou_ukb_nia_allvar_meta_analysis_chrposrefalt_cols_p_1e-5.TBL",show_col_types = F) %>%
    mutate(CHR = as.numeric(CHR))
  ukb_hits = vroom("../Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_p_1e-5.regenie",show_col_types = F) %>%
    mutate(`#CHROM` = as.numeric(`#CHROM`))
  aou_hits = vroom("../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos_gwsig.txt",show_col_types = F) %>%
    mutate(CHROM = as.numeric(CHROM))
  nia_hits = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos_gwsig.txt",show_col_types = F) %>%
    mutate(chr = as.numeric(chr))
  bed_file = data.frame(CHR=meta_hits$CHR,POSS=meta_hits$POS,POSE=meta_hits$POS) %>%
    add_row(data.frame(CHR=ukb_hits$`#CHROM`,POSS=ukb_hits$GENPOS,POSE=ukb_hits$GENPOS)) %>%
    add_row(data.frame(CHR=aou_hits$CHROM,POSS=aou_hits$GENPOS,POSE=aou_hits$GENPOS)) %>%
    add_row(data.frame(CHR=nia_hits$chr,POSS=nia_hits$pos,POSE=nia_hits$pos))
  vroom_write(bed_file,"working/all_study_hits.bed")
  print("Wrote BED file for AoU testing to: working/all_study_hits.bed")
}
{/R}

###########################
# BASH, on AoU
mkdir hwe_testing ; mkdir hwe_testing/hardy_out/
awk '{print $1 "\t" $2}' pgen_qc/chr1_geno_mac.psam > all_ids.txt
ancestries=(all eur afr amr sas eas mid)

# Run the hardy calculations
for ((chr=1;chr<=22;chr++)); do \
  for anc in "${ancestries[@]}"; do \
    ./plink2 --pfile pgen_qc/chr${chr}_geno_mac --keep ${anc}_ids.txt \
      --hardy midp --out hwe_testing/hardy_out/${chr}_${anc} --extract bed1 hwe_testing/all_study_hits.bed ;\
  done \
done

# Then merge the output
for anc in "${ancestries[@]}"; do \
  head -n 1 hwe_testing/hardy_out/1_${anc}.hardy > hwe_testing/hardy_out/${anc}_hwe_stats.txt ;\
  for file in hwe_testing/hardy_out/*${anc}.hardy; do \
    tail -n +2 "$file" >> hwe_testing/hardy_out/${anc}_hwe_stats.txt ;\
  done \
done

#############
# RCode on AoU to merge output into one table
all_var = vroom("hwe_testing/hardy_out/all_hwe_stats.txt",show_col_types = F) %>% 
    separate(ID,sep = "-",into=c("CHR",'POS'),remove = F) %>%
    mutate(IDrev = glue("{CHR}-{POS}-{AX}-{A1}")) %>% select(ID,IDrev,MIDP)
anc_hwe = data.frame(ID=all_var$ID,IDrev=all_var$IDrev,MIDP_ALL=NA,MIDP_AFR=NA,MIDP_AMR=NA,
                    MIDP_EAS=NA,MIDP_EUR=NA,MIDP_MID=NA,MIDP_SAS=NA)
hwe_datasets = lapply(X = c('all','afr','amr','eas','eur','mid','sas'),FUN = function(x) {
    vroom(glue("hwe_testing/hardy_out/{x}_hwe_stats.txt"),show_col_types = F)
})
for (row in 1:nrow(anc_hwe)) {
    curr_id = anc_hwe$ID[[row]] ; curr_id_rev = anc_hwe$IDrev[[row]]
    for (anc in c(1:7)) {
        curr_anc_id = which(hwe_datasets[[anc]]$ID == curr_id | hwe_datasets[[anc]]$ID == curr_id_rev)
        if (length(curr_anc_id) > 1) curr_anc_id = which(hwe_datasets[[anc]]$ID == curr_id)
        if (length(curr_anc_id) > 0)
            anc_hwe[[row,(2+anc)]] = hwe_datasets[[anc]]$MIDP[[curr_anc_id]]
    }
}
head(anc_hwe)
vroom_write(anc_hwe,"hwe_testing/anc_hwe_midp.txt")

print(nrow(anc_hwe))
print(nrow(anc_hwe %>% filter(MIDP_EUR <= 1e-15 | MIDP_AFR <= 1e-15 | MIDP_AMR <= 1e-15 |
                             MIDP_EAS <= 1e-15 | MIDP_SAS <= 1e-15 | MIDP_MID <= 1e-15)))
print(nrow(anc_hwe %>% filter(MIDP_EUR <= 1e-15 | MIDP_AFR <= 1e-15 | MIDP_AMR <= 1e-15)))

####################################
##############
# BACK TO FASRC






# To then give this file columns with single-ancestry regenie p-values, use the following:
ancestries=(all eur afr amr) ;\
for anc in "${ancestries[@]}"; do \
  echo $anc ;\
  awk 'NR==1 {$16 = "Pval"; $17 = "CHRPOS"; $18 = "IDrev"; print} NR>1 {$3 = $1 "-" $2 "-" $4 "-" $5; \
    $16 = 10^(-1 * $12); $17 = $1 "-" $2; $18 = $1 "-" $2 "-" $5 "-" $4; print}' \
    /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT/aou_AD_any_anc_${anc}_gwas.txt > \
    /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT/aou_AD_any_anc_${anc}_gwas_hyphenid_pval_chrpos.txt ;\
done

# Get data for multiancestry too
awk 'NR==1 {$18 = "IDrev"; print} NR>1 {$18 = $1 "-" $2 "-" $5 "-" $4; print}' \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos.txt > \
  /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT/aou_AD_any_anc_all_gwas_hyphenid_pval_chrpos.txt

# Intersect, checking for standard or reversed ID. First intersect is forwards
mkdir single_anc_final_hits_intersect ;\
for anc in "${ancestries[@]}"; do \
  echo $anc ;\
  if [ "$anc" == "all" ]; then
    awk 'NR==FNR{arr[$1];next} ($3 in arr) || ($17 in arr)' \
      /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/anc_hwe_midp.txt \
      /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT/aou_AD_any_anc_${anc}_gwas_hyphenid_pval_chrpos.txt > \
      single_anc_final_hits_intersect/final_hits_intersect_aou_${anc}.txt ;\
  else \
    awk 'NR==FNR{arr[$1];next} ($3 in arr) || ($18 in arr)' \
      /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/anc_hwe_midp.txt \
      /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT/aou_AD_any_anc_${anc}_gwas_hyphenid_pval_chrpos.txt > \
      single_anc_final_hits_intersect/final_hits_intersect_aou_${anc}.txt ;\
  fi ;\
done

# Then stick the p values into the table in R
{R}
make_supplemental_table_hwe_aou_p = function() {
  midp_df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/anc_hwe_midp.txt",show_col_types=F) %>%
    mutate(ALL_GW_P=NA,EUR_GW_P = NA,AFR_GW_P = NA,AMR_GW_P = NA)
  intersect_dfs = lapply(X=c("all","eur","afr","amr"),FUN = function(x) {
    vroom(glue("single_anc_final_hits_intersect/final_hits_intersect_aou_{x}.txt"),show_col_types = F)
  })
  for (r in 1:nrow(midp_df)) {
    if (r %% 100 == 0) print(glue("On row {r} of {nrow(midp_df)}"))
    curr_r = midp_df[r,]
    for (a in 1:4) {
      match = intersect_dfs[[a]] %>% filter(ID == curr_r$ID | IDrev == curr_r$ID)
      if (nrow(match) > 1) match %<>% filter(ID == curr_r$ID)
      if (nrow(match)>0) midp_df[[r,(8+a)]] = match$Pval
    }
  }
  write_xlsx(list("SupplementalTable" = midp_df), # remove empty rows
             path = "Paper_Tables_Figures/SupplementalTable_AoU_HWE_P_By_Ancestry.xlsx")
  print("Wrote table to: Paper_Tables_Figures/SupplementalTable_AoU_HWE_P_By_Ancestry.xlsx")
  
  return(midp_df)
}
{/R}
