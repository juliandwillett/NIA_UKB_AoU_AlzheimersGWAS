get_aou_only_fp = function() {
  # of all the GW significant hits in AoU, only 5 had a P val that became non sig with MCC
  # of these 5, only one was GW significant in the meta analysis (19-44899791-CTTTTTT-C)
  
  aou_hits = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/aou_allvar_anc_all_only_gw_sig.txt") %>%
    mutate(MCC_P=NA,MCC_Nonsig = F,ID = glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"))
  mcc_results = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/aou_AD_any_anc_all_gwas_mcc.txt") %>%
    mutate(ID = str_replace_all(ID,":","-"),ID = str_replace_all(ID,",","-")) %>%
    mutate(P=10^(-1 * LOG10P))
  
  # check p vals
  for (row in 1:nrow(aou_hits)) {
    print(row)
    mcc_match = mcc_results %>% filter(ID == aou_hits$ID[[row]])
    aou_hits$MCC_P[[row]] = mcc_match$P
    aou_hits$MCC_Nonsig[[row]] = ifelse(mcc_match$P > 5e-8,T,F)
  }
  vroom_write(aou_hits,"/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/aou_allvar_anc_all_only_gw_sig_with_mcc_vals.txt")
  return(aou_hits)
}
get_ukb_only_fp = function() {
  # of all the GW significant hits in AoU, only 5 had a P val that became non sig with MCC
  # of these 5, only one was GW significant in the meta analysis (19-44899791-CTTTTTT-C)
  
  ukb_hits = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab_withids_gwsig.txt") %>%
    mutate(MCC_P=NA,MCC_Nonsig = F,ID = glue("{`#CHROM`}-{GENPOS}-{ALLELE0}-{ALLELE1}")) %>%
    filter(`#CHROM` < 23)
  mcc_results = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/MCC_Results_1e-5/ukb_mcc_anc_all_merged.txt") %>%
    mutate(ID = glue("{`CHROM`}-{GENPOS}-{ALLELE0}-{ALLELE1}")) %>%
    mutate(P=10^(-1 * LOG10P))
  
  # check p vals
  for (row in 1:nrow(ukb_hits)) {
    print(row)
    mcc_match = mcc_results %>% filter(ID == ukb_hits$ID[[row]])
    ukb_hits$MCC_P[[row]] = mcc_match$P
    ukb_hits$MCC_Nonsig[[row]] = ifelse(mcc_match$P > 5e-8,T,F)
  }
  vroom_write(ukb_hits,"/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/Supplemental_Table_2_ukb_allvar_anc_all_only_gw_sig_with_mcc_vals.txt")
  return(aou_hits)
}
get_aou_pvals = function(df) {
  out_df = df
  
  # get data
  aou_summ_overlap = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/aou_AD_any_anc_all_gwas_tsv.txt") %>%
    mutate(CHRPOS = glue("{CHROM}-{GENPOS}")) %>% filter(CHRPOS %in% df$CHRPOS) %>%
    mutate(P=10^(-1 * LOG10P))
  
  # match alleles. Use for loop because not too big of a df
  for (row in 1:nrow(out_df)) {
    print(glue("On row {row}"))
    aou_match = aou_summ_overlap %>% filter(CHRPOS == out_df$CHRPOS[[row]])
    if (nrow(aou_match) > 0) { # check alleles
      aou_match2 = aou_match %>% filter((ALLELE0 == out_df$Allele0[[row]] & ALLELE1 == out_df$Allele1[[row]]) |
                                          (ALLELE0 == out_df$Allele1[[row]] & ALLELE1 == out_df$Allele0[[row]]))
      aou_match2_exact = aou_match %>% filter(ALLELE0 == out_df$Allele0[[row]] & ALLELE1 == out_df$Allele1[[row]])
      
      if (nrow(aou_match2) > 1) {
        out_df$AoU_Pval[[row]] = aou_match2_exact$P
        if (aou_match2_exact$BETA > 0) out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"+")
        else out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"-")
      }else if (nrow(aou_match2) > 0) {
        out_df$AoU_Pval[[row]] = aou_match2$P
        if (aou_match2$BETA > 0) out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"+")
        else out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"-")
      }else
        out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"?")
    }
  }
  return(out_df)
}
get_ukb_pvals = function(df) {
  out_df = df
  
  # get data
  ukb_summ_overlap = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/proxy_files_Step1_2_corrected_tab.txt") %>%
    mutate(CHRPOS = glue("{`#CHROM`}-{GENPOS}")) %>% filter(CHRPOS %in% out_df$CHRPOS) %>%
    mutate(P=10^(-1 * LOG10P))
  
  # match alleles. Use for loop because not too big of a df
  for (row in 1:nrow(out_df)) {
    print(glue("On row {row}"))
    ukb_match = ukb_summ_overlap %>% filter(CHRPOS == out_df$CHRPOS[[row]])
    if (nrow(ukb_match) > 0) { # check alleles
      ukb_match2 = ukb_match %>% filter((ALLELE0 == out_df$Allele0[[row]] & ALLELE1 == out_df$Allele1[[row]]) |
                                          (ALLELE0 == out_df$Allele1[[row]] & ALLELE1 == out_df$Allele0[[row]]))
      ukb_match2_exact = ukb_match %>% filter(ALLELE0 == out_df$Allele0[[row]] & ALLELE1 == out_df$Allele1[[row]])
      
      if (nrow(ukb_match2) > 1) {
        out_df$UKB_Pval[[row]] = ukb_match2_exact$P
        if (ukb_match2_exact$BETA > 0) out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"+")
        else out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"-")
      }else if (nrow(ukb_match2) > 0) {
        out_df$UKB_Pval[[row]] = ukb_match2$P
        if (ukb_match2$BETA > 0) out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"+")
        else out_df$Direction[[row]] = paste0(out_df$Direction[[row]],"-")
      }
    }
  }
  return(out_df)
}
get_aou_mcc = function(df) {
  out_df = df
  aou_mcc = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/aou_AD_any_anc_all_gwas_mcc.txt") %>%
    mutate(ID = str_replace_all(ID,":","-"),ID = str_replace_all(ID,",","-")) %>%
    mutate(P=10^(-1 * LOG10P))
  for (row in 1:nrow(out_df)) {
    mc_match = aou_mcc %>% filter(ID == out_df$ID[[row]])
    if (nrow(mc_match) > 0) out_df$AoU_MCC_P[[row]] = mc_match$P
  }
  return(out.df)
}
get_favor_annot = function(df,max_row) {
  out_df = data.frame(Locus=NA,ID=df$MarkerName,Rsid=NA,ProximalGene=NA,CHROM=df$CHR,POS=df$POS,
                      Allele0=toupper(df$Allele2),Allele1=toupper(df$Allele1),A1Freq=df$Freq1,
                      Beta=df$Effect,SE=df$StdErr,MetaP=df$`P-value`,Log10P=10^(-1*df$`P-value`),
                      Direction=NA,AoU_Pval=NA,UKB_Pval=NA,HetPVal=df$HetPVal,NewOld=NA,LZ_PostProb=NA,
                      CADD_Phred=NA,EnhancerRole=NA,EpiActivMax=NA,EpiReprMax=NA,
                      EpiTransMax=NA,Conserved=NA,FAVORAnnot=NA)
  out_df %<>% mutate(OR=round(exp(Beta),2),CI=glue("({round(exp(Beta-1.96*SE),2)}, {round(exp(Beta+1.96*SE),2)})"),.after="SE") %>%
    mutate(ID_FAVOR = glue("{CHROM}-{POS}-{Allele0}-{Allele1}")) %>%
    mutate(Direction = ifelse(Beta > 0,"+","-"))
  
  fuma_loci = vroom("locus_coords_grch38.bed",col_names=F) %>%
    rename(CHR=X1,START=X2,END=X3) %>% mutate(CHR = str_replace(CHR,"chr","")) %>%
    mutate(GenomicLocus = row_number()) %>% filter(CHR != "X")
  fuma_loci[c(14:46),'GenomicLocus'] = 13
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 46,GenomicLocus - 33, GenomicLocus))
  fuma_loci[c(90),'GenomicLocus'] = 56
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 90,GenomicLocus - 1, GenomicLocus))
  fuma_loci[c(93:94),'GenomicLocus'] = 58
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 94,GenomicLocus - 3, GenomicLocus))
  fuma_loci[c(135:137),'GenomicLocus'] = 98
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 137,GenomicLocus - 3, GenomicLocus))
  fuma_loci[c(150),'GenomicLocus'] = 110
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 150,GenomicLocus - 1, GenomicLocus))
  fuma_loci[c(156:164),'GenomicLocus'] = 115
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 164,GenomicLocus - 9, GenomicLocus))
  fuma_loci[c(172),'GenomicLocus'] = 122
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 172,GenomicLocus - 1, GenomicLocus))
  fuma_loci[c(173:175),'GenomicLocus'] = 123
  fuma_loci %<>% mutate(GenomicLocus = ifelse(row_number() > 175,GenomicLocus - 2, GenomicLocus))
  fuma_loci[c(29:38,40,42,44),'CHR'] = "2"
  fuma_loci[c(135:136),'CHR'] = "15"

  favor_annot = vroom("favor_hits_annotated_filtered_final.txt")
  bellenguez_loci = vroom("bellenguez_loci.txt")
  
  for (row in 1:nrow(out_df)) {
    print(glue("On row {row} of {nrow(out_df)}"))
    favor_match = which(favor_annot$`Variant (VCF)` == out_df$ID[[row]])
    fuma_locus_match = fuma_loci %>% filter(CHR == out_df$CHROM[[row]],START <= out_df$POS[[row]],
                                            END >= out_df$POS[[row]])
    old_loci = c(52,60:62,77,79,86,92:94,113,133,140:142,145:147,166,168)
    
    if (nrow(fuma_locus_match) > 0) {
      out_df$Locus[[row]] = fuma_locus_match$GenomicLocus
      if (fuma_locus_match$GenomicLocus %in% old_loci) out_df$NewOld[[row]] = "Old"
      else out_df$NewOld[[row]] = "New"
    }else{ # check against Bellenguez et al 2022
      curr_b_loc = bellenguez_loci %>% filter(CHR == out_df$CHROM[[row]],
                                              POS - 500000 <= out_df$POS[[row]],
                                              POS + 500000 >= out_df$POS[[row]])
      if (nrow(curr_b_loc) > 0) out_df$NewOld[[row]] = "Old"
      else out_df$NewOld[[row]] = "New"
    }
    
    if (!is.na(max_row) & row > max_row) break 
    
    if (length(favor_match) < 1) next # if not in favor database
    favor_row = favor_annot[favor_match,]
    print(glue("{favor_row$`Variant (VCF)`}"))
    rsid = favor_row$rsID
    prox_gene = glue("{favor_row$`Genecode Comprehensive Info`}")#; {favor_row$`UCSC Info`}")
    
    if (!is.na(favor_row$`CADD phred`)) out_df$CADD_Phred[[row]] = favor_row$`CADD phred`
    if (!is.na(favor_row$`CAGE Enhancer`) | !is.na(favor_row$GeneHancer) | !is.na(favor_row$SuperEnhancer))
      out_df$EnhancerRole[[row]] = "Yes"
    out_df$EpiActivMax[[row]] = max(c(favor_row$`aPC-Epigenetics-Active`,favor_row$H3K27ac,
                                      favor_row$H3K4me1,favor_row$H3K4me2,favor_row$H3K4me3,
                                      favor_row$H3K9ac,favor_row$H4k20me1,favor_row$H2AFZ),na.rm=T)
    out_df$EpiReprMax[[row]] = max(c(favor_row$`aPC-Epigenetics-Repressed`,favor_row$H3K9me3,
                                     favor_row$H3K27me3),na.rm=T)
    out_df$EpiTransMax[[row]] = max(c(favor_row$`aPC-Epigenetics-Transcription`,
                                      favor_row$H3K36me3,favor_row$H3K79me2),na.rm=T)
    
    # Full Favor annotation
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
    if (!is.na(favor_row$priPhCons) & favor_row$priPhCons >= 0.3) favor_sig %<>% paste(glue("Conservation priPhCons {favor_row$priPhCons};"))
    if (!is.na(favor_row$mamPhCons) & favor_row$mamPhCons >= 0.3) favor_sig %<>% paste(glue("Conservation mamPhCons {favor_row$mamPhCons};"))
    if (!is.na(favor_row$verPhCons) & favor_row$verPhCons >= 0.3) favor_sig %<>% paste(glue("Conservation verPhCons {favor_row$verPhCons};"))
    if (!is.na(favor_row$GerpN) & favor_row$GerpN >= 10) favor_sig %<>% paste(glue("Conservation GerpN {favor_row$GerpN};"))
    if (!is.na(favor_row$GerpS) & favor_row$GerpS >= 10) favor_sig %<>% paste(glue("Conservation GerpS {favor_row$GerpS};"))
    
    out_df$Rsid[[row]] = rsid
    out_df$ProximalGene[[row]] = prox_gene
    out_df$FAVORAnnot[[row]] = favor_sig
    if (str_detect(favor_sig,"Conservation")) out_df$Conserved[[row]] = "Yes"
    else out_df$Conserved[[row]] = "No"
  }
  return(out_df %>% mutate(CHRPOS = glue("{CHROM}-{POS}")))
}
getRsidsAndGenesForMissingVariants = function(df) { # obtained from gencode
  out_df = df %>% mutate(ProximalGene = str_replace_all(ProximalGene,"NONE\\(dist=NONE\\)",""))
  
  gencode = vroom("gencode.v44.annotation.gtf.gz",skip = 5,col_names = F) %>%
    filter(X3 == "gene") # focus on transcript to focus on TSS
  vep = vroom("missing_from_favor_vep.txt") %>% filter(str_detect(Existing_variation,'rs'))
  
  print("Started annotation")
  for (row in 1:nrow(out_df)) {
    if (row %% 500 == 0) print(glue("On row {row} of {nrow(out_df)}"))
    if (!is.na(out_df$ProximalGene[[row]])) next
    curr_row = out_df[row,]
    
    match_gencode = gencode %>% filter(X1 == paste0("chr",curr_row$CHROM),
                                       X4 >= curr_row$POS - 1e6, X4 <= curr_row$POS + 1e6) %>%
      separate(X9,into=c("GENE_ID","GENE_TYPE","GENE_NAME"),sep=";") %>%
      mutate(GENE_NAME = str_replace(GENE_NAME,'gene_name "',""),
             GENE_NAME = str_replace(GENE_NAME,'"',"")) %>%
      mutate(GENE_NAME = ifelse(X4 > curr_row$POS,glue("{GENE_NAME}(dist={X4-curr_row$POS})"),
                                ifelse(X4 < curr_row$POS,glue("{GENE_NAME}(dist={curr_row$POS-X4})"),
                                       GENE_NAME))) %>%
      mutate(Mindist = min(X4-curr_row$POS,curr_row$POS-X4)) %>%
      filter(!str_detect(GENE_NAME,"ENSG")) 
    if (nrow(match_gencode %>% filter(X4 <= curr_row$POS, X5 >= curr_row$POS)) > 0)
      match_gencode %<>% filter(X4 <= curr_row$POS, X5 >= curr_row$POS)
    else
      match_gencode %<>% arrange(Mindist) %>% slice(1:2)
    
    match_vep = vep %>% filter(str_detect(Location,glue("{curr_row$CHROM}:{curr_row$POS}")))
    
    if (nrow(match_gencode) < 1) out_df$ProximalGene[[row]] = "Noncoding"
    else out_df$ProximalGene[[row]] = paste(match_gencode$GENE_NAME,collapse=",")
    
    if (nrow(match_vep) > 0) out_df$Rsid[[row]] = match_vep$Existing_variation[[1]]
  }
  gc()
  return(out_df)
}
make_gsea_list = function(df) {
  genes = df %>% select(ProximalGene) %>% separate_rows(ProximalGene, sep = ",") %>%
    mutate(ProximalGene = str_replace_all(ProximalGene, "\\(.*?\\)", ""),
           ProximalGene = str_replace_all(ProximalGene," ",""))
  vroom_write(genes,"genes_for_web_gestalt.txt",col_names = F)
  print("Made file: genes_for_web_gestalt.txt")
}
