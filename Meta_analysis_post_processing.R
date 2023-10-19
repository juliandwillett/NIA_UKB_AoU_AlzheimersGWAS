# Filter meta analysis data for variants that demonstrate heterogeneity between studies
# Then write the file in a format that is readily plugged into LocusZoom

library(vroom)
library(tidyverse)

# Run this bash command in the terminal, system is too slow
# awk '{split($1, values, ":"); $(NF+1) = values[1]; $(NF+2) = values[2]; print $0}' \
aou_ukb_commonvar_meta_analysis.TBL > aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL

data = vroom("aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL") %>%
select(-`...17`) %>% rename(CHR = `MarkerName...16`, POS = `...18`)
data_qc_sorted = data %>% filter(HetPVal > 0.05) %>% arrange(CHR,POS)
vroom_write(data_qc_sorted,"aou_ukb_commonvar_meta_analysis_het_qc_sorted.txt")
