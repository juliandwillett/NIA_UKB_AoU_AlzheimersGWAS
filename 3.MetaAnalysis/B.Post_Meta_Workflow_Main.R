# Filter meta-analysis data for variants that demonstrate heterogeneity between studies
# Then write the file in a format that is readily plugged into LocusZoom

library(vroom)
library(tidyverse)
library(writexl)
library(glue)
library(magrittr)
library(readxl)
setwd("/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU")
`%notin%` <- function(x, y) !(x %in% y)
source('Functions.R')

# Run this bash command in the terminal, system is too slow
# awk '{split($1, values, ":"); $(NF+1) = values[1]; $(NF+2) = values[2]; print $0}' \
# aou_ukb_commonvar_meta_analysis.TBL > aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL

# this step requires 160 GB of RAM.
data = vroom("aou_ukb_allvar_meta_analysis_chrposrefalt.TBL") %>%
  select(-`...17`,-HetISq,-HetChiSq) %>% rename(CHR = `MarkerName...16`, POS = `...18`) %>% 
  filter(!(Effect >= 0 & str_detect(Direction,"\\-")),!(Effect < 0 & str_detect(Direction,"\\+"))) %>%
  filter(HetPVal > 5e-8, (MaxFreq - MinFreq) < 0.4) %>% arrange(CHR,POS) #QC per Bellenguez et al 2022
vroom_write(data,"aou_ukb_allvar_meta_het_freq_qc_sorted.txt") # Locus zoom file

# then assoc at GW significance level if GW significant in individual cohorts and meta analysis.

############################
# Now isolate GW significant hits for annotation
#gw_sig = data %>% filter(`P-value` <= 5e-8) %>% 
  #mutate(MarkerName...1 = str_replace_all(MarkerName...1,":","-"),MarkerName...1 = str_replace(MarkerName...1,",","-")) # so easily piped into favor
gw_sig = vroom("aou_ukb_allvar_meta_het_freq_qc_sorted.txt") %>% filter(`P-value` <= 5e-8) %>% 
  mutate(MarkerName = str_replace_all(MarkerName,":","-"),MarkerName = str_replace(MarkerName,",","-")) # so easily piped into favor
gw_sig %<>% filter(abs(Effect) < 5, CHR < 23)
vroom_write(gw_sig,"aou_ukb_allvar_meta_qc_sorted_gwsig.txt",col_names=T)

# Produce file for FAVOR batch annotator. 
vroom_write(gw_sig %>% select(MarkerName...1),"favor_hits.txt",col_names=F)

# Remove biased results
gw_sig_nomcc = vroom("aou_ukb_allvar_meta_qc_sorted_gwsig.txt") %>%
  filter(abs(Effect) < 5) %>% filter(MarkerName %notin% aou_hits_neg_mcc$ID)
vroom_write(gw_sig,"aou_ukb_allvar_meta_qc_sorted_gwsig_no_mccflags.txt")

############################
# Now start to annotate GW hits that are cleaned
gw_sig_nomcc = vroom("aou_ukb_allvar_meta_qc_sorted_gwsig_no_mccflags.txt")
df_hits = get_favor_annot(gw_sig_nomcc,NA) 
df_hits_aou_p = get_aou_pvals(df_hits) ; vroom_write(df_hits_aou_p,"tmp_meta_df.txt") ; gc()
df_hits_aou_ukb_p = get_ukb_pvals(df_hits_aou_p) ; vroom_write(df_hits_aou_ukb_p,"tmp_meta_df2.txt") ; gc()

# get info for variants missing from favor
# df_hits_all_annotated = getRsidsAndGenesForMissingVariants(df_hits_aou_ukb_p) ; vroom_write(df_hits_all_annotated,'meta_hits_all_annotated.txt')
df_hits_all_annotated = vroom("meta_hits_all_annotated.txt")
df_hits_proper_direction = df_hits_all_annotated %>% 
  filter(!(str_detect(Direction,"\\+") & str_detect(Direction,"\\-"))) %>%
  filter(!(is.na(AoU_Pval) & is.na(UKB_Pval))) 
df_hits_cadd = df_hits_proper_direction %>% filter(str_detect(FAVORAnnot,"CADD Phred"))
favor_annot = vroom("favor_hits_annotated_filtered_final.txt")
phred15 = df_hits_proper_direction %>% filter(ID %in% (favor_annot %>% filter(`CADD phred` >= 15))$`Variant (VCF)`) 
phred20 = df_hits_proper_direction %>% filter(ID %in% (favor_annot %>% filter(`CADD phred` >= 20))$`Variant (VCF)`) 

# make list for webgestalt and GSEA
make_gsea_list(df_hits_proper_direction)

###############
# Excel output

write_xlsx(df_hits_proper_direction, path = "meta_analysis_all_hits.xlsx")
write_xlsx(df_hits_cadd, path = "meta_analysis_cadd_hits.xlsx")
write_xlsx(phred15, path = "Table1-2_meta_analysis_cadd_hits_15.xlsx")
write_xlsx(phred20, path = "Table1-2_meta_analysis_cadd_hits_20.xlsx")


###############
# Get number of favor annot where at least one epigenetic field is greater than 10
nrow(favor_annot %>% filter(`Variant (VCF)` %in% df_hits_proper_direction$ID) %>%
       filter(H3K4me1 > 10 | H3K4me2 > 10 | H3K4me3 > 10 | H3K9ac > 10 | H3K9me3 > 10 | 
                H3K27ac > 10 | H3K27me3 > 10 | H3K36me3 > 10 | H3k79me2 > 10 | 
                H4k20me1 > 10 | H2AFZ > 10 | `aPC-Epigenetics-Active` > 10 | 
                `aPC-Epigenetics-Repressed` > 10 | `aPC-Epigenetics-Transcription` > 10 | 
                `aPC-Transcription-Factor` > 10))
