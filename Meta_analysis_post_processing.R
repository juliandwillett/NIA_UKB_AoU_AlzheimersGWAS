# Filter meta-analysis data for variants that demonstrate heterogeneity between studies
# Then write the file in a format that is readily plugged into LocusZoom

library(vroom)
library(tidyverse)
library(writexl)
library(glue)
library(magrittr)

# Run this bash command in the terminal, system is too slow
# awk '{split($1, values, ":"); $(NF+1) = values[1]; $(NF+2) = values[2]; print $0}' \
aou_ukb_commonvar_meta_analysis.TBL > aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL

data = vroom("aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL") %>%
select(-`...17`) %>% rename(CHR = `MarkerName...16`, POS = `...18`)
data_qc_sorted = data %>% filter(HetPVal > 0.05) %>% arrange(CHR,POS)
vroom_write(data_qc_sorted,"aou_ukb_commonvar_meta_analysis_het_qc_sorted.txt") # Locus zoom file

############################
# Now isolate GW significant hits for annotation
gw_sig = data_qc_sorted %>% filter(`P-value` <= 5e-8) %>% 
mutate(MarkerName...1 = str_replace_all(MarkerName...1,":","-"),MarkerName...1 = str_replace(MarkerName...1,",","-")) # so easily piped into favor

# Produce file for FAVOR batch annotator. 
vroom_write(gw_sig %>% select(MarkerName...1),"favor_hits.txt",col_names=F)

###########################
# Now take the FAVOR annotations, and make a file that readily works with the output in a well-organized file
favor_annot = vroom("aou_ukb_meta_favor_annotations.csv") %>% arrange(Chromosome,Position)

# Produce Excel file for writing annotations
df = data.frame(Locus=NA,ID=gw_sig$MarkerName...1,Rsid=NA,ProximalGene=NA,CHROM=gw_sig$CHR,POS=gw_sig$POS,
                Allele0=toupper(gw_sig$Allele2),Allele1=toupper(gw_sig$Allele1),A1Freq=gw_sig$Freq1,
                Beta=gw_sig$Effect,SE=gw_sig$StdErr,P=gw_sig$`P-value`,Log10P=10^(-gw_sig$`P-value`),
                HetPVal=gw_sig$HetPVal,NewOld=NA,FAVORAnnot=NA)
for (row in 1:nrow(df)) {
  print(glue("On row {row} of {nrow(df)}"))
  favor_match = which(favor_annot$`Variant (VCF)` == df$ID[[row]])
  if (length(favor_match) < 1) next # if not in favor database
  
  favor_row = favor_annot[favor_match,]
  rsid = favor_row$rsID
  prox_gene = glue("{favor_row$`Genecode Comprehensive Info`}; {favor_row$`UCSC Info`}")
  favor_sig = ""
  if (!is.na(favor_row$SuperEnhancer)) favor_sig %<>% append("SuperEnhancer; ")
  if (favor_row$H3K27me3 > 10) favor_sig %<>% append(glue("Transcription H3K27me3 {favor_row$H3K27me3}"))
  if (favor_row$H3K36me3 > 10) favor_sig %<>% append(glue("Transcription H3K36me3 {favor_row$H3K36me3}"))
  if (favor_row$totalRNA > 0.45) favor_sig %<>% append(glue("Transcription totalRNA {favor_row$totalRNA}"))
  if (row == 3) break
}
# will also want to replace UCSC transcript names with legible gene names

write_xlsx(df, path = "meta_analysis_hits.xlsx")
