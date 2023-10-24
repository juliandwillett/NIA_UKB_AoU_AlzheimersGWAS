# Filter meta-analysis data for variants that demonstrate heterogeneity between studies
# Then write the file in a format that is readily plugged into LocusZoom

library(vroom)
library(tidyverse)
library(writexl)
library(glue)
library(magrittr)

# Run this bash command in the terminal, system is too slow
# awk '{split($1, values, ":"); $(NF+1) = values[1]; $(NF+2) = values[2]; print $0}' \
# aou_ukb_commonvar_meta_analysis.TBL > aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL

#data = vroom("aou_ukb_commonvar_meta_analysis_IDcolon_chrposrefalt_cols.TBL") %>%
#  select(-`...17`) %>% rename(CHR = `MarkerName...16`, POS = `...18`)
#data_qc_sorted = data %>% filter(HetPVal > 0.05) %>% arrange(CHR,POS)
#vroom_write(data_qc_sorted,"aou_ukb_commonvar_meta_analysis_het_qc_sorted.txt") # Locus zoom file

############################
# Now isolate GW significant hits for annotation
#gw_sig = data_qc_sorted %>% filter(`P-value` <= 5e-8) %>% 
  #mutate(MarkerName...1 = str_replace_all(MarkerName...1,":","-"),MarkerName...1 = str_replace(MarkerName...1,",","-")) # so easily piped into favor
gw_sig = vroom("aou_ukb_commonvar_meta_analysis_het_qc_sorted.txt") %>% filter(`P-value` <= 5e-8) %>% 
  mutate(MarkerName = str_replace_all(MarkerName,":","-"),MarkerName = str_replace(MarkerName,",","-")) # so easily piped into favor


# Produce file for FAVOR batch annotator. 
vroom_write(gw_sig %>% select(MarkerName...1),"favor_hits.txt",col_names=F)

# Now take the FAVOR annotations, and make a file that readily works with the output in a well-organized file
favor_annot = vroom("aou_ukb_meta_favor_annotations.csv") %>% arrange(Chromosome,Position)

# Produce Excel file for writing annotations

get_favor_annot = function(max_row) {
  df = data.frame(Locus=NA,ID=gw_sig$MarkerName,Rsid=NA,ProximalGene=NA,CHROM=gw_sig$CHR,POS=gw_sig$POS,
                  Allele0=toupper(gw_sig$Allele2),Allele1=toupper(gw_sig$Allele1),A1Freq=gw_sig$Freq1,
                  Beta=gw_sig$Effect,SE=gw_sig$StdErr,MetaP=gw_sig$`P-value`,Log10P=10^(-gw_sig$`P-value`),
                  AoU_Pval=NA,UKB_Pval=NA,HetPVal=gw_sig$HetPVal,NewOld=NA,LZ_PostProb=NA,FAVORAnnot=NA)
  df %<>% mutate(OR=round(exp(Beta),2),CI=glue("({round(exp(Beta-1.96*SE),2)}, {round(exp(Beta+1.96*SE),2)})"),.after="SE") %>%
    mutate(ID_FAVOR = glue("{CHROM}-{POS}-{Allele0}-{Allele1}"))
  
  for (row in 1:nrow(df)) {
    print(glue("On row {row} of {nrow(df)}"))
    favor_match = which(favor_annot$`Variant (VCF)` == df$ID[[row]])
    if (length(favor_match) < 1) next # if not in favor database
    
    favor_row = favor_annot[favor_match,]
    
    print(glue("{favor_row$`Variant (VCF)`}"))
    if (!is.na(max_row) & row > max_row) break 
    
    rsid = favor_row$rsID
    prox_gene = glue("{favor_row$`Genecode Comprehensive Info`}")#; {favor_row$`UCSC Info`}")
    favor_sig = ""
    if (!is.na(favor_row$`CAGE Enhancer`)) favor_sig %<>% paste("CAGE Enhancer;")
    if (!is.na(favor_row$GeneHancer)) favor_sig %<>% paste(glue("Genehancer;"))
    if (!is.na(favor_row$SuperEnhancer)) favor_sig %<>% paste("SuperEnhancer;")
    if (!is.na(favor_row$`aPC-Conservation`) & favor_row$`aPC-Conservation` >= 10) favor_sig %<>% paste(glue("aPC-Conservation {favor_row$`aPC-Conservation`};"))
    if (!is.na(favor_row$`aPC-Epigenetics-Active`) & favor_row$`aPC-Epigenetics-Active` >= 10) favor_sig %<>% paste(glue("aPC-Epigenetics-Active {favor_row$`aPC-Epigenetics-Active`};"))
    if (!is.na(favor_row$`aPC-Epigenetics-Transcription`) & favor_row$`aPC-Epigenetics-Transcription` >= 10) favor_sig %<>% paste(glue("aPC-Epigenetics-Transcription {favor_row$`aPC-Epigenetics-Transcription`};"))
    if (!is.na(favor_row$`aPC-Epigenetics-Repressed`) & favor_row$`aPC-Epigenetics-Repressed` >= 10) favor_sig %<>% paste(glue("aPC-Epigenetics-Repressed {favor_row$`aPC-Epigenetics-Repressed`};"))
    if (!is.na(favor_row$`aPC-Local-Nucleotide-Diversity`) & favor_row$`aPC-Local-Nucleotide-Diversity` >= 10) favor_sig %<>% paste(glue("aPC-Local-Nucleotide-Diversity {favor_row$`aPC-Local-Nucleotide-Diversity`};"))
    if (!is.na(favor_row$`aPC-Transcription-Factor`) & favor_row$`aPC-Transcription-Factor` >= 10) favor_sig %<>% paste(glue("aPC-Transcription-Factor {favor_row$`aPC-Transcription-Factor`};"))
    if (!is.na(favor_row$`aPC-Mappability`) & favor_row$`aPC-Mappability` >= 10) favor_sig %<>% paste(glue("aPC-Mappability {favor_row$`aPC-Mappability`};"))
    if (!is.na(favor_row$`aPC-Mutation-Density`) & favor_row$`aPC-Mutation-Density` >= 10) favor_sig %<>% paste(glue("aPC-Mutation-Density {favor_row$`aPC-Mutation-Density`};"))
    if (!is.na(favor_row$`CADD phred`) & favor_row$`CADD phred` >= 10) favor_sig %<>% paste(glue("CADD Phred {favor_row$`CADD phred`};"))
    if (!is.na(favor_row$DNase) & favor_row$DNase >= 0.44) favor_sig %<>% paste(glue("Active DNase {favor_row$DNase};"))
    if (!is.na(favor_row$H2AFZ) & favor_row$H2AFZ >= 3.28) favor_sig %<>% paste(glue("Active H2AFZ {favor_row$H2AFZ};"))
    if (!is.na(favor_row$H3K4me1) & favor_row$H3K4me1 >= 4.5) favor_sig %<>% paste(glue("Active H3K4me1 {favor_row$H3K4me1};"))
    if (!is.na(favor_row$H3K4me2) & favor_row$H3K4me2 >= 3.5) favor_sig %<>% paste(glue("Active H3K4me2 {favor_row$H3K4me2};"))
    if (!is.na(favor_row$H3K4me3) & favor_row$H3K4me3 >= 3.7) favor_sig %<>% paste(glue("Active H3K4me3 {favor_row$H3K4me3};"))
    if (!is.na(favor_row$H3K9ac) & favor_row$H3K9ac >= 3.1) favor_sig %<>% paste(glue("Active H3K9ac {favor_row$H3K9ac};"))
    if (!is.na(favor_row$H3K9me3) & favor_row$H3K9me3 >= 3.7) favor_sig %<>% paste(glue("Repressed H3K9me3 {favor_row$H3K9me3};"))
    if (!is.na(favor_row$H3K27ac) & favor_row$H3K27ac >= 4.5) favor_sig %<>% paste(glue("Active H3K27ac {favor_row$H3K27ac};"))
    if (!is.na(favor_row$H3K27me3) & favor_row$H3K27me3 >= 4.5) favor_sig %<>% paste(glue("Repressed H3K27me3 {favor_row$H3K27me3};"))
    if (!is.na(favor_row$H3K36me3) & favor_row$H3K36me3 >= 3.9) favor_sig %<>% paste(glue("Transcription H3K36me3 {favor_row$H3K36me3};"))
    if (!is.na(favor_row$H3k79me2) & favor_row$H3k79me2 >= 3.5) favor_sig %<>% paste(glue("Transcription H3k79me2 {favor_row$H3k79me2};"))
    if (!is.na(favor_row$H4k20me1) & favor_row$H4k20me1 >= 3.7) favor_sig %<>% paste(glue("Active H4k20me1 {favor_row$H4k20me1};"))
    if (!is.na(favor_row$totalRNA) & favor_row$totalRNA >= 0.1) favor_sig %<>% paste(glue("Transcription totalRNA {favor_row$totalRNA};"))
    if (!is.na(favor_row$priPhyloP) & favor_row$priPhyloP >= 0.3) favor_sig %<>% paste(glue("Conservation priPhyloP {favor_row$priPhyloP};"))
    if (!is.na(favor_row$priPhCons) & favor_row$priPhCons >= 0.9) favor_sig %<>% paste(glue("Conservation priPhCons {favor_row$priPhCons};"))
    if (!is.na(favor_row$mamPhCons) & favor_row$mamPhCons >= 0.9) favor_sig %<>% paste(glue("Conservation mamPhCons {favor_row$mamPhCons};"))
    if (!is.na(favor_row$verPhCons) & favor_row$verPhCons >= 0.9) favor_sig %<>% paste(glue("Conservation verPhCons {favor_row$verPhCons};"))
    if (!is.na(favor_row$GerpN) & favor_row$GerpN >= 10) favor_sig %<>% paste(glue("Conservation GerpN {favor_row$GerpN};"))
    if (!is.na(favor_row$GerpS) & favor_row$GerpS >= 10) favor_sig %<>% paste(glue("Conservation GerpS {favor_row$GerpS};"))
        
    df$Rsid[[row]] = rsid
    df$ProximalGene[[row]] = prox_gene
    df$FAVORAnnot[[row]] = favor_sig
  }
  
  return(df)
}

df_hits = get_favor_annot(NA)
df_hits_cadd = df_hits %>% filter(str_detect(FAVORAnnot,"CADD Phred"))

write_xlsx(df_hits, path = "meta_analysis_all_hits.xlsx")
write_xlsx(df_hits_cadd, path = "meta_analysis_cadd_hits.xlsx")
