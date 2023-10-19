# Filter meta analysis data for variants that demonstrate heterogeneity between studies
# Then write the file in a format that is readily plugged into LocusZoom

library(vroom)
library(tidyverse)
library(writexl)

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

# Produce file for FAVOR batch annotator
vroom_write(gw_sig %>% select(MarkerName...1),"favor_hits.txt",col_names=F)

# Produce Excel file for writing annotations
df = data.frame(Locus=NA,ID=gw_sig$MarkerName...1,Rsid=NA,ProximalGene=NA,CHROM=gw_sig$CHR,POS=gw_sig$POS,
                Allele0=toupper(gw_sig$Allele2),Allele1=toupper(gw_sig$Allele1),A1Freq=gw_sig$Freq1,
                Beta=gw_sig$Effect,SE=gw_sig$StdErr,P=gw_sig$`P-value`,Log10P=10^(-gw_sig$`P-value`),
                HetPVal=gw_sig$HetPVal,NewOld=NA,FAVORAnnot=NA)
write_xlsx(df, path = "meta_analysis_hits.xlsx")
