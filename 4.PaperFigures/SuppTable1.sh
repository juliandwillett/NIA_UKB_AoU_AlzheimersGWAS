# Code for obtaining numbers for Supplemental Table 1

# SUPP TABLE 1: INFO
# NIAGADS
retained_nia = vroom("/n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/keep_selfALL25660.txt",delim = '\t',col_names = F)
nia_dem = vroom("../Data_Links/NIAGADS_Personal/NIAGADS_demographic_data_all.txt") %>%
  filter(IID %in% retained_nia$X1)
nia_dem_defined_ages = nia_dem %>% filter(!is.na(Age))
table(nia_dem$Affection.Status)
table(nia_dem$Affection.Status)/25660
table(nia_dem$superpop2)
table(nia_dem$superpop2)/25660
table((nia_dem %>% filter(Affection.Status == 1))$superpop2)
c(mean(nia_dem_defined_ages$Age),sd(nia_dem_defined_ages$Age))
table(nia_dem$Sex)/25660
gw_sig = vroom("../Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_p_1e-5.txt") %>%
  filter(chr == 19,pos == 44908684 | pos == 44908822)

# UKB
ukb_pheno_covar = vroom("../Data_Links/UKB_GWAS_Data/regenie_pheno_covar.txt")
table(ukb_pheno_covar$ad_proxy)
table(ukb_pheno_covar$ad_proxy)/159629
ukb_anc = vroom("../Data_Links/UKB_GWAS_Data/ancestry_group_stats")
ukb_anc
c(mean(ukb_pheno_covar$Age_at_recruitment),sd(ukb_pheno_covar$Age_at_recruitment))
table(ukb_pheno_covar$Sex)/159629
gw_sig = vroom("../Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_gwsig.regenie") %>%
  filter(`#CHROM` == 19,GENPOS == 44908684 | GENPOS == 44908822)

# AOU: RUN ON WORKBENCH
aou_covar = vroom("regenie_input/regenie_covar_commonpcs_with_anc.txt",show_col_types = F)
aou_pheno = vroom("regenie_input/regenie_pheno.txt",show_col_types = F)
merged = merge(aou_covar,aou_pheno,by='IID')
table(merged$AD)
table(merged$AD_any)
nrow(merged)
table(merged$AD_any) / 244838
table(merged$ancestry_pred)
table(merged$ancestry_pred) / 244838
table((merged %>% filter(AD_any == 1))$ancestry_pred)
c(mean(merged$Age),sd(merged$Age))
table(merged$Sex) / 244838

# FOR VARIANT AF, COLLECT ON FASRC
gw_sig = vroom("../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_p_1e-1.txt") %>%
  filter(`CHROM` == 19,GENPOS == 44908684 | GENPOS == 44908822)
