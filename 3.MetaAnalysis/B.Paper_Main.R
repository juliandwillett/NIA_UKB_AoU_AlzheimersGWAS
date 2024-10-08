library(vroom)
library(tidyverse)
library(writexl)
library(glue)
library(magrittr)
library(readxl)
library(ggrepel)
library(viridis)
setwd("/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU")
`%notin%` <- function(x, y) !(x %in% y)
source('Rcode/Functions.R')
options(warn = 1)

########################
### ORGANIZE SINGLE STUDY DATA

# GET VARIANTS/BEDFILES FOR HWE THEN INTERSECT ALL THE DATASETS
# bed_coords = make_files_for_hwe()
# intersect_ids = get_nominal_by_gw_sig(get_for_intersect = T)

# ORGANIZE BASE DATA
# df_all_hits = get_nominal_by_gw_sig(get_for_intersect = F,start_index = 1)
# vroom_write(df_all_hits,"working/df_all_hits_across_studies.txt")

# ORGANIZE DATA FOR CLUMPING ON AOU (MORE DIVERSE)
# prep_data_for_separate_loci_allcohorts(df_all_hits)

# ASSIGN HITS TO LOCI, GET FAVOR GENE ANNOT
# df_nia_hits = get_favor_annot(df_all_hits %>% filter(NIA_P<=5e-8))
# df_nimh_hits = get_favor_annot(df_all_hits %>% filter(NIMH_P<=5e-8))
# df_nia_nimh_meta_hits = get_favor_annot(df_all_hits %>% filter(NIA_NIMH_META_P<=5e-8))
# df_ukb_hits = get_favor_annot(df_all_hits %>% filter(UKB_P<=5e-8))
# df_aou_hits = get_favor_annot(df_all_hits %>% filter(AOU_P<=5e-8))
# df_ukb_aou_meta_hits = get_favor_annot(df_all_hits %>% filter(UKB_AOU_META_P<=5e-8))
# gwas_hits_annotated = list(df_nia_hits,df_nimh_hits,df_nia_nimh_meta_hits,df_ukb_hits,
#                            df_aou_hits,df_ukb_aou_meta_hits)
# names(gwas_hits_annotated) = c("NIA","NIMH","NIA_NIMH","UKB","AOU","UKB_AOU")
# saveRDS(gwas_hits_annotated,"working/single_gwas_annot.rds")

# SEPARATE LOCI USING CLUMPING
df_all_hits = vroom("working/df_all_hits_across_studies.txt")
gwas_hits_annotated = readRDS("working/single_gwas_annot.rds")
df_hits_ind_loci = separate_loci_gwas(df_all_hits,gwas_hits_annotated) # clumped loci for all datasets
df_hits_ind_loci_clin_meta_nominal = get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=T,
                                                            get_direction = T)

# # ADD CELL EXPR CHANGE ANNOTATIONS
# clin_ad_diffexpr = add_cell_expression(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META,p_cutoff = 0.01)
# clin_ad_diffexpr = addDiseaseGroup(clin_ad_diffexpr,p_cutoff = 0.01)
# 
# ad_by_proxy_diffexpr = add_cell_expression(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META,p_cutoff = 0.01)
# ad_by_proxy_diffexpr = addDiseaseGroup(ad_by_proxy_diffexpr,p_cutoff = 0.01)
# 
# saveRDS(list(clin_ad_diffexpr,ad_by_proxy_diffexpr),"working/ad_adbyproxy_diffexpr.rds")

qc_expr_data = readRDS("working/ad_adbyproxy_diffexpr.rds")
names(qc_expr_data) = c("AD","AD_by_proxy")

###################################
# COUNTING FOR PAPER
# NUM VARIANTS/SAMPLES ORIGINALLY IN NIAGADS
# wc -l /n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/NIAGADS_pgen/cc34438_combined.pvar
# wc -l /n/holystore01/LABS/tanzi_lab/Users/dmitry/NIAGADS_v9_analysis/NIAGADS_pgen/cc34438_combined.psam

# NUM VARIANTS NIAGADS POST QC, IE IN FINAL DATA SUMM STATS. RUN IN TERMINAL
# awk '$1 != "X" && NR > 1 {print}' /n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_hwe_chrpos.txt | wc -l
max(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_p_gwsig.txt")$N,na.rm = T)

# NUM VARIANTS/SAMPLES/CASES UKB POST QC
# awk 'NR>1 && $1 != 23 {print}' /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete.regenie | wc -l
max(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_p_1e-5.regenie")$N,na.rm = T)
table(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/regenie_pheno_covar.txt")$ad_proxy)

# NUM VARIANTS/SAMPLES AOU POST QC
# awk 'NR>1 {print}' /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs.txt | wc -l
max(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_p_1e-5.txt")$N)

# NUMBER OF LOCI IN UKB PRE CHECKING FOR NOMINAL
ukb_pre = get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,
                                             get_direction = T)$UKB
tmp = getNumberCommonRareLoci(ukb_pre,maf_cut = 0.01,return_lead=T,cohort = "UKB")

# NUMBER OF LOCI IN AOU PRE CHECKING FOR NOMINAL
aou_pre = get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,
                                 get_direction = T)$AOU
tmp = getNumberCommonRareLoci(aou_pre,maf_cut = 0.01,return_lead=T,cohort = "AOU")

# NUMBER OF GENES IN AD BY PROXY META FINAL SET
count_genes(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META %>% match_signs_proxy() %>% remove_multiallelics())

# NUMBER OF LOCI IN UKB AOU META ONLY DETECTED IN ONE DATASET
only_aou = count_reproduced_loci(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META %>% match_signs_proxy() %>% remove_multiallelics(),c('UKB_P'))
ukb_and_aou = count_reproduced_loci(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META %>% remove_multiallelics(),c("UKB_P",'AOU_P'))

####################################
### FIGURES AND TABLES

####################################
### FIGURE 1: STUDY SUMMARY

# Num participants in studies
nia_summ = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_p_1e-5.txt")
max(nia_summ$N)
ukb_summ = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_p_1e-5.regenie")$N
aou_summ = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_p_1e-5.txt")$N
print(max(ukb_summ) + max(aou_summ))

# NIA SIG VARIANTS AND LOCUS COUNTS
tmp = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$NIA,maf_cut = 0.01,
                              return_lead=T,cohort = "NIA")

# NIA-NIMH SIG VARIANTS AND LOCUS COUNTS
tmp = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META,
                              maf_cut = 0.01,return_lead=T,cohort = "NIA_NIMH_META")

# UKB SIG VARIANTS AND LOCUS COUNTS
tmp = getNumberCommonRareLoci(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$UKB),
                              maf_cut = 0.01,return_lead=T,cohort = "UKB")

# AOU SIG VARIANTS AND LOCUS COUNTS
tmp = getNumberCommonRareLoci(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$AOU) %>%
                                remove_multiallelics(),
                              maf_cut = 0.01,return_lead=T,cohort = "AOU")

# UKB AOU META SIG VARIANTS AND LOCUS COUNTS
tmp = getNumberCommonRareLoci(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META),
                              maf_cut = 0.01,return_lead=T,cohort = "UKB_AOU_META")

# GET UKB/AOU SIGNAL OVERLAP WITH CLINICAL AD
# OLD CHECKING OF OVERLAP OF SIGNALS
# nrow(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META %>% filter(UKB_AOU_META_P <= 0.05)) / 
#   nrow(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META)
# tmp = get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,get_direction = T)
# ukb_stat = nrow(tmp$UKB %>% filter(ClinNom)) / nrow(tmp$UKB)
# ukb_stat
# aou_stat = nrow(tmp$AOU %>% filter(ClinNom)) / nrow(tmp$AOU)
# aou_stat
# ukb_aou_meta_stat = nrow(tmp$UKB_AOU_META %>% filter(ClinNom)) / nrow(tmp$UKB_AOU_META)
# ukb_aou_meta_stat
# CHECK BY LOCI, SO ANY VARIANT IN EACH LOCUS
nia_nimh_overlap = count_reproduced_loci(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META,"UKB_AOU_META_P")
length(unique(nia_nimh_overlap$Locus)) / 22
ukb_overlap = count_reproduced_loci(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,
                                             get_direction = T)$UKB,"NIA_NIMH_META_P") %>% match_signs_proxy()
length(unique(ukb_overlap$Locus)) / 30
aou_overlap = count_reproduced_loci(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,
                                             get_direction = T)$AOU,"NIA_NIMH_META_P") %>% match_signs_proxy() %>% remove_multiallelics()
length(unique(aou_overlap$Locus)) / 1558
ukb_aou_overlap = count_reproduced_loci(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,
                                             get_direction = T)$UKB_AOU_META,"NIA_NIMH_META_P") %>% match_signs_proxy() %>% remove_multiallelics()
length(unique(ukb_aou_overlap$Locus)) / 1345

# GET OVERLAP WITH FAVOR AND NUMBER THAT OVERLAP/LINK TO GENE ENHANCERS
nrow(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META))
nrow(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META %>% filter(EnhancerRole == "Yes")))

nrow(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META))
nrow(match_signs_proxy(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META %>% filter(EnhancerRole == "Yes")))

# GET NUMBER OF VARIANTS IN GENES THAT ARE DIFFERENTIALLY EXPRESSED
tmp = make_table_cell_expression(df = qc_expr_data$AD,col = 'NIA_NIMH_META',
                                 clean_cell_columns = T,include_signs=T,only_new=F,enhancer_role = F)

tmp = make_table_cell_expression(qc_expr_data$AD_by_proxy,col = "UKB_AOU_META",
                                 clean_cell_columns = T,include_signs=T,only_new=F,enhancer_role = F)

########################
# FIGURE 2: PRODUCED ON AOU

########################
# FIGURE 3: MANHATTANS GWAS AND META ANALYSIS OF CLINICAL AD
lead_nia_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$NIA,
                                       maf_cut = 0.01,return_lead=T,cohort = "NIA")
nia_man = make_manhattan_easy_label("NIAGADS",df_hits_ind_loci_clin_meta_nominal$NIA,
                                    lead_var = lead_nia_all,fig_num = "3A",hide_new_gene=F,hide_x=T)
lead_nia_nimh_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META,
                                            maf_cut = 0.01,return_lead=T,cohort = "NIA_NIMH_META")
nia_nimh_meta_man = make_manhattan_easy_label("NIA_NIMH_META",
                                              df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META,
                                    lead_var = lead_nia_nimh_all,fig_num = "3B",hide_new_gene=F,
                                    only_common = F,
                                    hide_x = T,no_indels = F)

#########################
# TABLE 1: GWAS LEAD VAR NIAGADS + META LEAD VAR NIAGADS+NIMH
lead_nia_nimh_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$NIA_NIMH_META,
                                            maf_cut = 0.01,return_lead=T,cohort = "NIA_NIMH_META")
table1 = make_table_lead_loci("NIA_NIMH",lead_nia_nimh_all)
write_xlsx(table1,"Paper_Tables_Figures/table1.xlsx")

##########################
# FIGURE 4: MANHATTAN AOU/UKB-AOU
# This is now S2B
lead_aou_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$AOU,
                                       maf_cut = 0.01,return_lead=T,cohort = "AOU")
lead_aou_all_matching_clinAD_dir = remove_multiallelics(match_signs_proxy(lead_aou_all))
aou_man = make_manhattan_easy_label("AOU",df_hits_ind_loci_clin_meta_nominal$AOU %>%
                                      remove_multiallelics() %>% match_signs_proxy(),
                                    lead_var = lead_aou_all_matching_clinAD_dir,fig_num = "4A_dirmatch",hide_new_gene=F,
                          only_common = F,hide_x = T,no_indels = T,no_multiallelics=T,
                          maf_cut = NA)

lead_ukb_aou_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META,
                                       maf_cut = 0.01,return_lead=T,cohort = "UKB_AOU_META") 
lead_ukb_aou_all_dir = getNumberCommonRareLoci(match_signs_proxy(lead_ukb_aou_all),
                                               maf_cut = 0.01,return_lead=T,cohort = "UKB_AOU_META")
ukb_aou_man = make_manhattan_easy_label(dataset="UKB_AOU",df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META %>%
                                          remove_multiallelics() %>% match_signs_proxy(),
                                    lead_var = lead_ukb_aou_all_dir,fig_num = "4B",hide_new_gene = F,
                          only_common = F,hide_x= T,no_indels = T,no_multiallelics = T)

###########################
# TABLE 2: LEAD VARIANTS UKB-AOU META. NOT PUTTING AOU ALONE RE GENERALIZABILITY
lead_ukb_aou_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$UKB_AOU_META,
                                           maf_cut = 0.01,return_lead=T,cohort = "UKB_AOU_META")
lead_ukb_aou_all_dir = getNumberCommonRareLoci(match_signs_proxy(lead_ukb_aou_all),
  maf_cut = 0.01,return_lead=T,cohort = "UKB_AOU_META")
table2 = make_table_lead_loci("UKB_AOU",lead_ukb_aou_all_dir) 
write_xlsx(table2,"Paper_Tables_Figures/table2.xlsx")

###########################
# TABLES 3-4: SINGLE CELL EXPRESSION OUTCOMES
# GENERAL OUTCOMES
table3a = make_table_cell_expression(df = qc_expr_data$AD,col = "NIA_NIMH_META",
                                    clean_cell_columns = T,include_signs=T,only_new=T,
                                    enhancer_role = F) #%>% filter(MAF >= 0.01)
table3b = match_signs_proxy(make_table_cell_expression(qc_expr_data$AD_by_proxy,col = "UKB_AOU_META",
                                    clean_cell_columns = T,include_signs=T,only_new=T,
                                    enhancer_role = F))
table3 = table3a %>% add_row() %>% add_row(table3b)

write_xlsx(table3,path = "Paper_Tables_Figures/table3.xlsx")

# PATH SPECIFIC OUTCOMES
table4a = make_table_cell_expression_bypath(qc_expr_data$AD,coh_col = "NIA_NIMH_META",clean_cell_columns = T,
                                                 include_signs=T,only_new=T,enhancer_linked = F) %>%
  filter(Path %in% c("NFT","PlaqN"))#,MAF >= 0.01)
table4b = make_table_cell_expression_bypath(qc_expr_data$AD_by_proxy,coh_col = "UKB_AOU_META",clean_cell_columns = T,
                                            include_signs=T,only_new=T,enhancer_linked = F) %>%
  filter(Path %in% c("NFT","PlaqN"))#,MAF >= 0.01)
table4 = table4a %>% add_row() %>% add_row(table4b)
write_xlsx(table4,path = "Paper_Tables_Figures/table4.xlsx")

#######################
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

# AOU - run on AOU workbench, given privacy policy, using similar code, except for getting MAFs
gw_sig = vroom("../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_p_1e-1.txt") %>%
  filter(`CHROM` == 19,GENPOS == 44908684 | GENPOS == 44908822)

#######################
# SUPP TABLE 2: SUMMARY STATS ACROSS STUDIES
st2 = make_df_all_hits_table(df_all_hits %>% filter(CHR != 23)) %>%
  mutate(CHR = as.numeric(CHR)) %>% arrange(CHR) %>%
  relocate(CHR,.after=AOU_EUR_DIR)
write_xlsx(st2,"Paper_Tables_Figures/supp_table2.xlsx")

########################
# SUPP TABLE 3: IND LOCI IN CLINICAL AD
st3 = make_supp_table_all_loci_variants(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=T,get_direction = T)$NIA_NIMH_META,
                                        'NIA_NIMH_META') %>%
  arrange(CHR,ID) %>% select(-Flipped) %>% rename(VariantID=ID,Chromosome=CHR,
                                                 EffectAllele=EA) %>%
  relocate(VariantID,.after="Chromosome") %>%
  relocate(Direction,.after = "NIMH_P") %>%
  make_consecutive_loci()
write_xlsx(st3,path = "Paper_Tables_Figures/supp_table3.xlsx")

########################
# SUPP TABLE 4: IND LOCI IN AD-by-proxy
st4 = make_supp_table_all_loci_variants(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,get_direction = T)$UKB_AOU_META,
                                        'UKB_AOU_META') %>%
  arrange(CHR,ID) %>% rename(Chromosome = CHR,VariantID = ID) %>% select(-Flipped) %>%
  make_consecutive_loci()
write_xlsx(st4,path = "Paper_Tables_Figures/supp_table4.xlsx")

########################
# SUPP TABLE 5-6: IND LOCI IN UKB AND AOU
st5 = make_supp_table_all_loci_variants(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,get_direction = T)$UKB,
                                        'UKB_AOU_META') %>%
  arrange(CHR,ID) %>% rename(Chromosome = CHR,VariantID = ID) %>% select(-Flipped) %>%
  make_consecutive_loci()
write_xlsx(st5,path = "Paper_Tables_Figures/supp_table5.xlsx")

st6 = make_supp_table_all_loci_variants(get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,get_direction = T)$AOU,
                                        'UKB_AOU_META') %>%
  arrange(CHR,ID) %>% rename(Chromosome = CHR,VariantID = ID) %>% select(-Flipped) %>%
  make_consecutive_loci()
write_xlsx(st6,path = "Paper_Tables_Figures/supp_table6.xlsx")

########################
# SUPP FIG 2: AOU MANHATTAN WO INDELS
lead_aou_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$AOU,
                                       maf_cut = 0.01,return_lead=T,cohort = "AOU") %>%
  filter(nchar(EA) == nchar(NEA),nchar(EA) == 1)
lead_aou_all_matching_clinAD_dir = match_signs_proxy(lead_aou_all)
aou_man = make_manhattan_easy_label("AOU",df_hits_ind_loci_clin_meta_nominal$AOU %>%
                                      filter(str_sub(UKB_AOU_DIRECTION,2,2) == str_sub(NIA_NIMH_DIRECTION,1,1)),
                                    lead_var = lead_aou_all_matching_clinAD_dir,fig_num = "S2",hide_new_gene=F,
                                    only_common = F,hide_x = T,no_indels=T,hide_all_genes=F,
                                    no_multiallelics = T)

########################
# SUPP FIGURE 3: MANHATTAN UKB
lead_ukb_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$UKB,
                                       maf_cut = 0.01,return_lead=T,cohort = "UKB")
lead_ukb_all_matching_clinAD_dir = match_signs_proxy(lead_ukb_all)
make_manhattan_easy_label("UKB",df_hits_ind_loci_clin_meta_nominal$UKB,
                          lead_var = lead_ukb_all_matching_clinAD_dir,fig_num = "S3",hide_new_gene=F,
                          only_common = F,hide_x = T,
                          no_indels = T,no_multiallelics = T)

########################
# SUPP TABLE 7: BELLENGUEZ CLIN AD VS AD BY PROXY 
# supp_table = comparePValsBellenguezLoci(make_for_testing = T)
st7 = comparePValsBellenguezLoci(make_for_testing = F) %>% select(-IDrev)
nrow(st7 %>% filter(AOU_P <= 0.05,ifelse(Bell_OR>=1,AOU_OR>=1,AOU_OR<=1)))
write_xlsx(st7,"Paper_Tables_Figures/supp_table7.xlsx")

######################
# SUPP FIGURE 4: MANHATTAN ANC SPECIFIC
lead_aou_all = getNumberCommonRareLoci(df_hits_ind_loci_clin_meta_nominal$AOU,
                                       maf_cut = 0.01,return_lead=T,cohort = "AOU")
lead_aou_all_matching_clinAD_dir = remove_multiallelics(match_signs_proxy(lead_aou_all))
make_manhattan_easy_label("AOU_AFR",df_hits_ind_loci_clin_meta_nominal$AOU %>%
                            filter(str_sub(UKB_AOU_DIRECTION,2,2) == str_sub(NIA_NIMH_DIRECTION,1,1)),
                          lead_var = lead_aou_all_matching_clinAD_dir,fig_num = "S4A",hide_new_gene=F,
                          only_common = F,hide_x = T,
                          no_indels = T,no_multiallelics = T)
make_manhattan_easy_label("AOU_AMR",df_hits_ind_loci_clin_meta_nominal$AOU %>%
                            filter(str_sub(UKB_AOU_DIRECTION,2,2) == str_sub(NIA_NIMH_DIRECTION,1,1)),
                          lead_var = lead_aou_all_matching_clinAD_dir,fig_num = "S4B",hide_new_gene=F,
                          only_common = F,hide_x = T,
                          no_indels = T,no_multiallelics = T)
make_manhattan_easy_label("AOU_EUR",df_hits_ind_loci_clin_meta_nominal$AOU %>%
                            filter(str_sub(UKB_AOU_DIRECTION,2,2) == str_sub(NIA_NIMH_DIRECTION,1,1)),
                          lead_var = lead_aou_all_matching_clinAD_dir,fig_num = "S4C",hide_new_gene=F,
                          only_common = F,hide_x = T,
                          no_indels = T,no_multiallelics = T)

########################
# OTHER STATISTICS
# Number of AoU loci when removing indels
multiallelics = vroom("working/aou_multiallelics.txt",delim="\t")
aou_no_indels_no_multiallelics = get_multi_nominal_loci(df_hits_ind_loci,get_clin_nom=F,get_direction = T)$AOU %>% 
  filter(nchar(EA) == nchar(NEA),nchar(EA) == 1) %>% mutate(CHRPOS = glue("{CHR}-{POS}")) %>%
  filter(CHRPOS %notin% multiallelics$CHRPOS)
tmp = getNumberCommonRareLoci(aou_no_indels_no_multiallelics,maf_cut = 0.01,
                              return_lead=T,cohort = "AOU")
(640)/(1558)
# to get count of loci in AoU that were nominal and not indel/multi, just eyeball 
#   lead_aou_all_matching_dir variable since there are only 10 rows

# Number of AoU loci when removing indels AND monovariant loci
aou_no_indels_no_monovariant_loci = aou_no_indels_no_multiallelics %>% group_by(Locus) %>%
  filter(n() > 1) %>% ungroup(Locus) 
tmp = getNumberCommonRareLoci(aou_no_indels_no_monovariant_loci,maf_cut = 0.01,
                              return_lead=T,cohort = "AOU")

# Organized STs
write_xlsx(list("ST2" = st2,"ST3" = st3,"ST4" = st4,"ST5" = st5,"ST6" = st6,
                "ST7" = st7),"Paper_Tables_Figures/supplemental_tables.xlsx")

