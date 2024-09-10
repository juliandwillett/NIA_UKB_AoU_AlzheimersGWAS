get_inflation = function(dataset) {
  case = NA ; control = NA
  if (dataset == "NIA") {
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_hwe_chrpos.txt") %>%
      filter(!is.na(p),chr != "X")
    case = 10565 ;  control = 25660 - case
    lambda = median(qchisq(df$p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
    lambda1000 = 1 + (lambda - 1) * (1 / case + 1 / control) * 500
  }else if (dataset == "UKB") {
    # awk 'NR==1 || $1 != 23 {print $11}' /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete.regenie > /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_justchisq_autosomes.regenie
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_justchisq_autosomes.regenie",delim='\t') 
    case = 25785 ; control = 159628 - case
    lambda = median(df$CHISQ,na.rm=T)/qchisq(0.5,1)
    lambda1000 = 1 + (lambda - 1) * (1 / case + 1 / control) * 500
  }else if (dataset == 'AOU') {
    # awk '{print $11}' /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs.txt > /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_justchisq_autosomes.txt
    df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_justchisq_autosomes.txt",delim='\t') 
    case = 11290 ; control = 244845 - case
    lambda = median(df$CHISQ,na.rm=T)/qchisq(0.5,1)
    lambda1000 = 1 + (lambda - 1) * (1 / case + 1 / control) * 500
  }else if (dataset == 'UKB_AOU') {
    df = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/final_meta_results/meta_ukb_aou_justpvals.TBL",delim='\t') %>%
      filter(!is.na(`P-value`),`P-value` >= 0,`P-value` <= 1)
    case = (25785+11290) ; control = (159628 + 244845) - case
    lambda = median(qchisq(df$`P-value`, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
    lambda1000 = 1 + (lambda - 1) * (1 / case + 1 / control) * 500
  }
  
  lambda = median(df$CHISQ,na.rm=T)/qchisq(0.5,1)
  lambda = median(qchisq(p, df=1, lower.tail=FALSE)) / qchisq(0.5, 1)
  lambda1000 = 1 + (lambda - 1) * (1 / case + 1 / control) * 500
}
getUKBDataforFigure2 = function() {
  ukb_samples = vroom("ukb_psam_wgs.psam")
  ukb_covar = vroom("~/true_lab_storage/Data_Links/UKB_Covar_PhenoFiles/covar2") %>% 
    filter(IID %in% ukb_samples$IID)
  ukb_pheno = vroom("~/true_lab_storage/Data_Links/UKB_Covar_PhenoFiles/pheno2") %>% 
    filter(IID %in% ukb_samples$IID)
  ukb_anc = vroom("~/true_lab_storage/Data_Links/UKB_Covar_PhenoFiles/keep_white.tsv") %>% 
    filter(IID %in% ukb_samples$IID)
  
  table_for_aou = (merge(ukb_covar,ukb_pheno,by="IID") %>% mutate(White=(IID %in% ukb_anc$IID)) %>%
    select(-IID,-FID.x,-FID.y))
  vroom_write(table_for_aou[-c(3:23),],"ukb_covar_pheno_for_fig_2.txt")
}
get_aou_pvals = function(df) {
  out_df = df
  
  # get data
  aou_summ_overlap = vroom("working/meta_hits_aou_intersect.txt",col_names=F)
  names(aou_summ_overlap) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos_gwsig.txt"))
  names(aou_summ_overlap)[15:16] = c("P","CHRPOS")

  # match alleles. Use for loop because not too big of a df
  for (row in 1:nrow(out_df)) {
    if (row %% 250 == 0) print(glue("On row {row}"))
    r = out_df[row,]
    aou_match = aou_summ_overlap %>% mutate(ID = str_replace_all(ID,":","-"),ID = str_replace(ID,",","-")) %>%
      filter(ID == r$ID) 
    if (nrow(aou_match) > 0) { out_df$AoU_Pval[[row]] = aou_match$P } 
  }
  # compare stats to flipped allele (which varies)
  tmp = out_df %>% filter(is.na(AoU_Pval)) %>% mutate(Flipped_ID = glue("{CHROM}-{POS}-{Allele1}-{Allele0}"))
  flipped_present = which(tmp$Flipped_ID %in% aou_summ_overlap$ID)
  aou_pval_flipped = aou_summ_overlap$P[which(aou_summ_overlap$ID %in% tmp$Flipped_ID)]
  print(glue("Flipped p-vals: Meta {tmp$MetaP[flipped_present]}. AoU {aou_pval_flipped}"))
  
  return(out_df)
}
get_ukb_pvals = function(df) {
  # flipped p to use by row number to deal with flipped alleles
  out_df = df
  
  # get data
  ukb_summ_overlap = vroom("working/meta_hits_ukb_intersect.txt",col_names=F) 
  names(ukb_summ_overlap) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_gwsig.regenie"))
  names(ukb_summ_overlap)[[15]] = "CHRPOS"

  # match alleles. Use for loop because not too big of a df
  for (row in 1:nrow(out_df)) {
    if (row %% 250 == 0) print(glue("On row {row}"))
    r = out_df[row,]
    ukb_match = ukb_summ_overlap %>% filter(ID == r$ID)
    if (nrow(ukb_match) > 1 & length(unique(ukb_match$Pval)) == 1) { out_df$UKB_Pval[[row]] = ukb_match$Pval[[1]] }
    else if (nrow(ukb_match) > 0) { out_df$UKB_Pval[[row]] = ukb_match$Pval }
  }
  
  # compare stats to flipped allele (which varies)
  tmp = out_df %>% filter(is.na(UKB_Pval)) %>% mutate(Flipped_ID = glue("{CHROM}-{POS}-{Allele1}-{Allele0}"))
  flipped_present = which(tmp$Flipped_ID %in% ukb_summ_overlap$ID)
  ukb_pval_flipped = ukb_summ_overlap$Pval[which(ukb_summ_overlap$ID %in% tmp$Flipped_ID)]
  print(glue("Flipped p-vals: Row_num {flipped_present}. ID {tmp$ID[flipped_present]}. Meta {tmp$MetaP[flipped_present]}. AoU {tmp$AoU_Pval[flipped_present]}. UKB Normal {}. UKB Flip {ukb_pval_flipped}"))
  
  return(out_df %>% mutate(ProximalGene = ifelse(ProximalGene == "NA",NA,ProximalGene),
                           ProximalGene = str_replace(ProximalGene,"NONE\\(dist=NONE\\)","")))
}
get_niagads_pvals = function(df) {
  out_df = df %>% mutate(Niagads_Pval = NA,.after="UKB_Pval")
  
  # get data
  nia_int = vroom("working/meta_hits_niagads_intersect.txt",col_names = F)
  names(nia_int) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/NIAGADS_meta_chrpos_gwsig.txt"))

  # match alleles. Use for loop because not too big of a df
  for (row in 1:nrow(out_df)) {
    if (row %% 250 == 0) print(glue("On row {row}"))
    r = out_df[row,]
    match = nia_int %>% mutate(marker = str_replace_all(marker,":","-"),marker = str_replace(marker,",","-")) %>%
      filter(marker == r$ID) 
    if (nrow(match) > 0) { out_df$Niagads_Pval[[row]] = match$p } 
  }
  # compare stats to flipped allele (which varies)
  tmp = out_df %>% filter(is.na(Niagads_Pval)) %>% mutate(Flipped_ID = glue("{CHROM}:{POS}:{Allele1}:{Allele0}"))
  flipped_present = which(tmp$Flipped_ID %in% nia_int$marker)
  pval_flipped = nia_int$p[which(nia_int$marker %in% tmp$Flipped_ID)]
  print(glue("Flipped p-vals: Meta {tmp$MetaP[flipped_present]}. Niagads {pval_flipped}"))
  
  return(out_df)
}
# get_bellenguez_pvals = function(df) {
#   out_df = df %>% mutate(Bellenguez_Pval = NA,.after="Niagads_Pval")
#   
#   # get data
#   bell_int = vroom("meta_hits_bell_intersect_aou_vs_ukb.txt",col_names = F)
#   names(bell_int) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/Other_GWAS/bellengeuz_meta_ids_chrpos_gwsig.tsv"))
#   
#   # match alleles. Use for loop because not too big of a df
#   for (row in 1:nrow(out_df)) {
#     if (row %% 250 == 0) print(glue("On row {row}"))
#     r = out_df[row,]
#     match = bell_int %>% mutate(ID = str_replace_all(ID,":","-"),ID = str_replace(ID,",","-")) %>%
#       filter(ID == r$ID) 
#     if (nrow(match) > 0) { out_df$Bellenguez_Pval[[row]] = match$Pval } 
#   }
#   # compare stats to flipped allele (which varies)
#   tmp = out_df %>% filter(is.na(Bellenguez_Pval)) %>% mutate(Flipped_ID = glue("{CHROM}:{POS}:{Allele1}:{Allele0}"))
#   flipped_present = which(tmp$Flipped_ID %in% bell_int$ID)
#   pval_flipped = bell_int$Pval[which(bell_int$ID %in% tmp$Flipped_ID)]
#   print(glue("Flipped p-vals: Meta {tmp$MetaP[flipped_present]}. Bellenguez {pval_flipped}"))
#   
#   return(out_df)
# }
getKnownLoci = function() {
  # # for getting TSS of known loci from Rudy's list
  # gencode = vroom("gencode.v44.annotation.gtf.gz",skip = 5,col_names = F) %>%
  #   separate(X9,into=c("GENE_ID","GENE_TYPE","GENE_NAME"),sep=";") %>%
  #   mutate(GENE_NAME = str_replace(GENE_NAME,'gene_name "',""),
  #          GENE_NAME = str_replace(GENE_NAME,'"',"")) %>%
  #   mutate(GENE_NAME = str_replace(GENE_NAME," ","")) %>%
  #   mutate(X1 = str_replace(X1,"chr","")) %>%
  #   filter(X3 == "gene") # focus on transcript to focus on TSS, which is X4
  # all_loci = read_xlsx(list.files(".",pattern="CAF")) %>% drop_na(`Locus - Lead Gene (Bold)`) %>%
  #   rename(GENE_NAME = `Locus - Lead Gene (Bold)`) %>%
  #   mutate(GENE_NAME = ifelse(GENE_NAME == "MAGAT5","MGAT5",GENE_NAME),
  #          GENE_NAME = ifelse(GENE_NAME == "AC099552.4","ENSG00000217825",GENE_NAME),
  #          GENE_NAME = ifelse(GENE_NAME == "PAX1P","PAX1",GENE_NAME),
  #          GENE_NAME = ifelse(GENE_NAME == "MTSS1L","MTSS2",GENE_NAME),
  #          GENE_NAME = ifelse(GENE_NAME == "BZRAP1-AS1","TSPOAP1-AS1",GENE_NAME),
  #          GENE_NAME = ifelse(GENE_NAME == "ATPB8","ATP8B1",GENE_NAME),
  #          GENE_NAME = ifelse(GENE_NAME == "C15ORF41","CDIN1",GENE_NAME)) %>%
  #   add_row(data.frame(GENE_NAME = gencode[which(str_detect(gencode$GENE_NAME,"MS4A")),'GENE_NAME'])) %>%
  #   add_row(data.frame(GENE_NAME = c("MAF","WWOX")))
  # known_internal = merge(all_loci,gencode,by="GENE_NAME") %>%
  #   add_row(data.frame(GENE_NAME = "AC087500.1",X1="17",X4=5192068)) %>% # not in gencode, select transcript
  #   add_row(data.frame(GENE_NAME = "AC004687.2",X1="17",X4=58331728)) %>% #wightman et al hit
  #   add_row(data.frame(GENE_NAME = "NTN6",X1="19",X4=49164664)) %>% # internal data
  #   mutate(CHRPOS = glue("{X1}-{X4}"))
  # for (g in 208:250) { # checking for imbalances (all are more than 500 kb from TSS)
  #   m = df_hits %>% filter(str_detect(ProximalGene,known_internal$GENE_NAME[[g]]),NewOld == "New")
  #   if (nrow(m)>0)
  #     print(glue("GeneTSS {known_internal$CHRPOS[[g]]}. Variant {m$CHRPOS}"))
  # }
  # 
  # bellenguez_data = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/Other_GWAS/bellengeuz_meta_ids.tsv") %>%
  #   filter(Pval <= 5e-8) %>% mutate(CHRPOS = glue("{chromosome}-{base_pair_location}"))
  # wightman_data = vroom("/n/holystore01/LABS/tanzi_lab/Users/jwillett/Data_Links/Other_GWAS/Wightmanetal2021_PGCALZ2ExcludingUKBand23andME_METALInverseVariance_MetaAnalysis.txt.gz") %>%
  #   filter(p_value <= 5e-8) %>% mutate(X4 = glue("{chromosome}:{base_pair_location+1}")) # GRCh37
  # # vroom_write(data.frame(chr=paste0("chr",wightman_data$chromosome),pos=wightman_data$base_pair_location,pos=wightman_data$base_pair_location + 1),
  # #                        "wightman_for_liftover_grch37.txt",col_names = F)
  # wightman_data_grch38_gwsig = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/Other_GWAS/wightman_for_liftover_grch38_gwsig.bed",col_names = F) %>%
  #   mutate(X4 = sub("-.*", "", X4)) %>% mutate(X4 = str_replace(X4,"chr",""))
  # wightman_grch38_data = merge(wightman_data,wightman_data_grch38_gwsig,by="X4") %>%
  #   mutate(CHRPOS = glue("{X1}-{X2}"),CHRPOS = str_replace(CHRPOS,"chr","")) %>% select(-X4,-base_pair_location,-chromosome)
  # 
  # # now merge the three and output it. Only need CHR and POS
  # out_df = data.frame(CHRPOS = known_internal$CHRPOS) %>%
  #   add_row(data.frame(CHRPOS = bellenguez_data$CHRPOS)) %>%
  #   add_row(data.frame(CHRPOS = wightman_grch38_data$CHRPOS)) %>%
  #   separate(CHRPOS,into = c("CHR","POS"),sep = "-")
  # 
  # vroom_write(out_df,"known_loci_chrpos.txt")
  
  return(vroom("known_loci_chrpos.txt"))
}
assignLociNumbers = function(df) {
  # Number the loci
  curr_loc = 1
  out_df = df
  for (row in 1:nrow(df)) {
    print(row)
    r = out_df[row,]
    if (row %% 1000 == 0) print(glue("On row {row} of {nrow(df)}"))
    if (row == 1) { out_df$Locus[[row]] = curr_loc ; next }
    r_prev = out_df[row-1,]
    if (is.na(r$Locus)) {
      if (row == 1) { print("Code must be revised. NA locus in first row") ; stop() } # for other analyses, not this one
      if (r$CHR == r_prev$CHR & (r$POS-r_prev$POS <= 500000)) { # joined with prior locus
        out_df$Locus[[row]] = r_prev$Locus
      }else{
        curr_loc = curr_loc + 1
        out_df$Locus[[row]] = curr_loc
      }
    }
  }
  print(glue("Num total loci pre-QC: {length(unique(out_df$Locus))}"))
  return(out_df)
}
get_favor_annot = function(df,window=500,maf_cut=0.01) {
  # Window is the NewOld locus definition in kb
  # Matches by EXACT ID MATCH, not ID or IDrev match
  
  # Makes df gets requisite data
  out_df = data.frame(Locus=NA,ID=df$ID,IDrev=df$IDrev,Rsid=NA,ProximalGene=NA,CHROM=df$CHR,POS=df$POS,
                      Allele0=toupper(df$NEA),Allele1=toupper(df$EA),NewOld=NA,
                      EnhancerRole=NA,EnhancerLinkedGene=NA,CADD_Phred=NA)
  favor_annot = vroom("../Data_Links/FAVOR/favor_hits_annot.txt",skip=1,delim=",",col_names=F)
  names(favor_annot) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",n_max=1,delim="\t"))
  favor_annot %<>% rename(ID = variant_vcf)
  known_loci = getKnownLoci() #GW sig or likely causal
  
  favor_rsids = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/FAVOR_Rsids.txt") # space delim
  super_enhancers = vroom("dbSUPER_SuperEnhancers_hg19.tsv")
  
  # Go row by row to get all required fields
  for (row in 1:nrow(out_df)) {
    print(glue("On row {row} of {nrow(out_df)}"))
    favor_row = match_accounting_for_indels(favor_annot,out_df[row,]) %>% distinct(.keep_all = T)

    curr_b_loc = known_loci %>% filter(CHR == out_df$CHROM[[row]],
                                            POS - window*1000 <= out_df$POS[[row]],
                                            POS + window*1000 >= out_df$POS[[row]]) 
    if (nrow(curr_b_loc) > 0) out_df$NewOld[[row]] = "Old"
    else out_df$NewOld[[row]] = "New"
    
    if (nrow(favor_row) < 1) {  # if not in favor database
      print(glue("{out_df$ID[[row]]} not in database"))
      if (length(which(favor_annot$ID == out_df$ID[[row]]))<1)
        if (!(str_detect(out_df$ID[[row]],"\\*")))
          stop("Missing entries from FAVOR. Need to rerun checking the database. Will place empty annotations for missing entries, so error does not recur.")
      next
    }
    prox_gene = glue("{favor_row$genecode_comprehensive_info[[1]]}")
    
    if (!is.na(favor_row$cadd_phred[[1]])) out_df$CADD_Phred[[row]] = favor_row$cadd_phred[[1]]
    if (!is.na(favor_row$cage_enhancer[[1]]) | !is.na(favor_row$genehancer[[1]]) | !is.na(favor_row$super_enhancer[[1]])) {
      out_df$EnhancerRole[[row]] = "Yes"
      if (!is.na(favor_row$genehancer[[1]])) {
        # The score is a factor, but what is more important is the existence of a genehancer that provides
        # evidence that there is some relationship between the annotation and a gene's expression. The
        # scores are also not necessarily interpretable since they can diverge by so much.
        lgenes = regmatches(favor_row$genehancer[[1]], gregexpr("(?<=connected_gene=)[^;]+",favor_row$genehancer[[1]], perl = TRUE))[[1]]
        out_df$EnhancerLinkedGene[[row]] = paste(lgenes,collapse=",")
        # print("Gene linked to variant by GeneHancer")
      }
      if (!is.na(favor_row$super_enhancer[[1]])) {
        senhancers = str_split(favor_row$super_enhancer[[1]],",")[[1]]
        db_match = which(super_enhancers$se_id %in% senhancers)
        gene_matches = unique(super_enhancers$gene_symbol[db_match])
        out_df$EnhancerLinkedGene[[row]] = paste(gene_matches,collapse=",")
      }
    }
    out_df$ProximalGene[[row]] = prox_gene

    if (out_df$ID[[row]] %in% favor_rsids$FAVOR_VCF) {
      out_df$Rsid[[row]] = unique(favor_rsids$Rsid[which(favor_rsids$FAVOR_VCF == out_df$ID[[row]])])
    }
  }
  
  # Number the loci
  out_df = assignLociNumbers(out_df %>% mutate(POS = as.numeric(POS)))
  
  return(out_df %>% mutate(CHRPOS = glue("{CHROM}-{POS}")))
}
get_rsids_only = function(df) {
  rsids = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/FAVOR_Rsids.txt",show_col_types = F)
  out_df = df
  for (row in 1:nrow(df)) {
    match = rsids %>% filter(FAVOR_VCF == df$ID[[row]])
    if (nrow(match)>0) out_df$Rsid[[row]] = match$Rsid[[1]]
  }
  return(out_df)
}
getRsidsAndGenesForMissingVariants = function(df,only_closest) { # obtained from gencode
  print("To execute bash script for getting the annotations from the large files using dbsnp, \\
        use this file: /n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/getRsids.sh. \\
        Preferable to use FAVOR batch iterator, though this is here.")
  
  out_df = df %>% mutate(ProximalGene = str_replace_all(Gene,"NONE\\(dist=NONE\\)",""))
  
  gencode = vroom("gencode.v44.annotation.gtf.gz",skip = 5,col_names = F,show_col_types = F) %>%
    filter(X3 == "gene") # focus on transcript to focus on TSS
  vep_gene = vroom("missing_from_favor_vep.txt",show_col_types = F)
  favor_rsids = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/FAVOR_Rsids.txt",show_col_types = F)
  vep_rsids = vroom("missing_from_favor_vep.txt",show_col_types = F) %>% 
    filter(str_detect(Existing_variation,"rs"))
  
  print("Started annotation")
  for (row in 1:nrow(out_df)) {
    if (row %% 500 == 0) print(glue("On row {row} of {nrow(out_df)}"))
    curr_row = out_df[row,]
    
    if (!str_detect(curr_row$ID,"\\*")) { # not spanning deletion
      fav = favor_rsids %>% filter(FAVOR_VCF == curr_row$ID_FAVOR) %>% distinct(.keep_all=T)
      if (nrow(fav) > 0) {out_df$Rsid[[row]] = unique(fav$Rsid)}
      else { # not in favor, check in VEP
        ve = vep_rsids %>% filter(str_detect(Location,glue("{curr_row$CHROM}:{curr_row$POS}")))
        if (nrow(ve) > 0) out_df$Rsid[[row]] = unique(ve$Existing_variation)
        else print(curr_row$ID)
      }
    }

    # Get Gene information
    if (!is.na(out_df$ProximalGene[[row]])) next
    match_gencode = gencode %>% filter(X1 == paste0("chr",curr_row$CHROM),
                                       X4 >= curr_row$POS - 1e6, X4 <= curr_row$POS + 1e6) %>%
      separate(X9,into=c("GENE_ID","GENE_TYPE","GENE_NAME"),sep=";") %>%
      mutate(GENE_NAME = str_replace(GENE_NAME,'gene_name "',""),
             GENE_NAME = str_replace(GENE_NAME,'"',"")) %>%
      mutate(GENE_NAME = ifelse(X4 > curr_row$POS,glue("{GENE_NAME}(dist={X4-curr_row$POS})"),
                                ifelse(X4 < curr_row$POS,glue("{GENE_NAME}(dist={curr_row$POS-X4})"),
                                       GENE_NAME))) %>%
      mutate(Mindist = X4-curr_row$POS,curr_row$POS-X4) %>%
      filter(!str_detect(GENE_NAME,"ENSG")) %>% arrange(X1,X4)
    if (nrow(match_gencode %>% filter(X4 <= curr_row$POS, X5 >= curr_row$POS)) > 0)
      match_gencode %<>% filter(X4 <= curr_row$POS, X5 >= curr_row$POS) 
    if (nrow(match_gencode) > 0) {
      if (only_closest)
        out_df$ProximalGene[[row]] = paste((match_gencode %>% slice_min(abs(Mindist), n = 1))$GENE_NAME,collapse=",")
      else
        out_df$ProximalGene[[row]] = paste((match_gencode)$GENE_NAME,collapse=",")
    }
  }
  gc()
  return(out_df %>% mutate(ProximalGene = ifelse(ProximalGene == "-()",NA,ProximalGene)) %>%
           mutate(ProximalGene = str_replace_all(ProximalGene,"()",""),
                  ProximalGene = str_replace_all(ProximalGene,"(dist=-)","")))
}
make_gsea_list = function(df,maf_cut = 0.01) {
  print(glue("Maf cut: {maf_cut}"))
  common = df %>% filter(A1Freq >= maf_cut)
  rare = df %>% filter(A1Freq < maf_cut)
  
  genes_c = common %>% select(ProximalGene) %>% add_row(ProximalGene = df$EnhancerLinkedGene) %>% 
    mutate(ProximalGene = str_replace_all(ProximalGene, "\\(dist.*?\\)", ""),
           ProximalGene = str_replace_all(ProximalGene," ",""),
           ProximalGene = gsub("\\([^)]*\\)", "", ProximalGene)) %>%
    separate_rows(ProximalGene, sep = ",") %>% distinct()

  vroom_write(genes_c,"genes_for_web_gestalt_common.txt",col_names = F)
  
  genes_r = rare %>% select(ProximalGene) %>% add_row(ProximalGene = df$EnhancerLinkedGene) %>% 
    mutate(ProximalGene = str_replace_all(ProximalGene, "\\(dist.*?\\)", ""),
           ProximalGene = str_replace_all(ProximalGene," ",""),
           ProximalGene = gsub("\\([^)]*\\)", "", ProximalGene)) %>%
    separate_rows(ProximalGene, sep = ",") %>% distinct()
  
  vroom_write(genes_r,"genes_for_web_gestalt_rare.txt",col_names = F)
  
  vroom_write(genes_c %>% add_row(genes_r) %>% distinct(),"genes_for_web_gestalt_all.txt",col_names=F)
  
  print("Made file: genes_for_web_gestalt.txt")
  return(list(genes_c,genes_r))
}
comparePValsBellenguezLoci = function(make_for_testing) {
  bell_hits_pval = read_xlsx("Bellenguez_2022_SuppTables.xlsx",sheet = "Supplementary Table 5",skip=2) %>%
    select(Varianta,Chr.b,Positionc,Gened,`Minor/Major allele`,MAFe,`...13`,`Stage I + II`) %>% slice(-1) %>%
    separate(`Minor/Major allele`, into = c("MINOR", "MAJOR"), sep = "/") %>%
    separate(`Stage I + II`,into = c("OR"),sep="\\(") %>%
    rename_all(~ c("Rsid","CHR","POS","GENE","MINOR","MAJOR","MAF","Pval_Final","OR")) %>%
    select(-GENE) %>%
    mutate(ID=glue("{CHR}-{POS}-{MINOR}-{MAJOR}")) %>% #POS GRCh38
    mutate(IDrev = glue("{CHR}-{POS}-{MAJOR}-{MINOR}")) 
  
  if (make_for_testing) {
    hwe_test = data.frame(CHR=bell_hits_pval$CHR,POSS=bell_hits_pval$POS,POSE=bell_hits_pval$POS)
    vroom_write(hwe_test,"working/bell_hits_for_hwe_test.bed",col_names = F)
    print("Wrote bed file to: working/bell_hits_for_hwe_test.bed")
    print("Run HWE stats using: https://github.com/juliandwillett/AoU_AlzheimersGWAS/blob/main/1.AoU_GWAS/K.CheckForFalsePositives_Anc.sh")
    return()
  }
  
  # vroom_write(data.frame(ID=bell_hits_pval$ID) %>% add_row(ID=bell_hits_pval$IDrev),
  #             "bellenguez_hits_id.txt")
  
  # Will intersect these hits with meta results to save time on computational time
  # INTERSECT CLINICAL AD HWE
  # awk 'NR==FNR{arr[$1]; next} $1 in arr' bellenguez_hits_id.txt \
  # final_meta_results/meta_niagads_nimh_hwe.TBL > \
  # working/bell_clinad_intersect_hwe.txt
  
  # INTERSECT CLINICAL AD NO-HWE
  # awk 'NR==FNR{arr[$1]; next} $1 in arr' bellenguez_hits_id.txt \
  # final_meta_results/meta_niagads_nimh_nohwe.TBL > \
  # working/bell_clinad_intersect_no_hwe.txt
  
  # INTERSECT UKB
  # awk 'NR==FNR{arr[$1]; next} $3 in arr' bellenguez_hits_id.txt \
  # /n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id.regenie > \
  # working/bell_ukb_intersect.txt
  
  # INTERSECT AOU
  # awk 'NR==FNR{arr[$1]; next} $3 in arr' bellenguez_hits_id.txt \
  # /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs.txt > \
  # working/bell_aou_intersect.txt
  
  # INTERSECT AOU AFR
  # awk 'NR==FNR{arr[$1]; next} $3 in arr' bellenguez_hits_id.txt \
  # /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_afr_summ_stats_AD_any_allchr.txt > \
  # working/bell_aou_afr_intersect.txt
  
  # INTERSECT AOU AMR
  # awk 'NR==FNR{arr[$1]; next} $3 in arr' bellenguez_hits_id.txt \
  # /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_amr_summ_stats_AD_any_allchr.txt > \
  # working/bell_aou_amr_intersect.txt
  
  # INTERSECT AOU EUR
  # awk 'NR==FNR{arr[$1]; next} $3 in arr' bellenguez_hits_id.txt \
  # /n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_eur_summ_stats_AD_any_allchr.txt > \
  # working/bell_aou_eur_intersect.txt
  
  # INTERSECT AD BY PROXY META
  # awk 'NR==FNR{arr[$1]; next} $1 in arr' bellenguez_hits_id.txt \
  # final_meta_results/meta_ukb_aou.TBL > \
  # working/bell_adbyproxy_intersect.txt
  
  out_df = data.frame(ID = bell_hits_pval$ID,CHR=bell_hits_pval$CHR,Rsid=bell_hits_pval$Rsid,
                      EA=bell_hits_pval$MINOR,Bell_OR=bell_hits_pval$OR,
                      Bell_P=bell_hits_pval$Pval_Final,AD_P=NA,UKB_P=NA,AOU_OR=NA,
                      AOU_P=NA,AOU_AFR_P=NA,AOU_AMR_P=NA,AOU_EUR_P=NA,AD_by_proxy_P=NA,HWE_ALL=NA,
                      HWE_AFR=NA,HWE_AMR=NA,HWE_EAS=NA,HWE_EUR=NA,HWE_MID=NA,
                      HWE_SAS=NA,IDrev=bell_hits_pval$IDrev)
  ad_data = vroom("working/bell_clinad_intersect_no_hwe.txt",show_col_types = F,col_names = F) %>%
    rename(ID = X1,P = X10)
  ukb_data = vroom("working/bell_ukb_intersect.txt",show_col_types = F,col_names = T) %>%
    mutate(P = 10^(-LOG10P))
  aou_data = vroom("working/bell_aou_intersect.txt",show_col_types = F,col_names = T) %>%
    mutate(P = 10^(-LOG10P), MAF = A1FREQ,Effect = BETA,Direction = sign(Effect)) 
  aou_data = make_all_minor_alleles(aou_data,direction_index = c(1))
  
  # get data for single ancestries
  aou_data_afr = vroom("working/bell_aou_afr_intersect.txt",show_col_types = F,col_names = T) %>%
    mutate(P = 10^(-LOG10P), MAF = A1FREQ,Effect = BETA,Direction = sign(Effect)) 
  aou_data_afr = make_all_minor_alleles(aou_data_afr,direction_index = c(1))
  aou_data_amr = vroom("working/bell_aou_amr_intersect.txt",show_col_types = F,col_names = T) %>%
    mutate(P = 10^(-LOG10P), MAF = A1FREQ,Effect = BETA,Direction = sign(Effect)) 
  aou_data_amr = make_all_minor_alleles(aou_data_amr,direction_index = c(1))
  aou_data_eur = vroom("working/bell_aou_eur_intersect.txt",show_col_types = F,col_names = T) %>%
    mutate(P = 10^(-LOG10P), MAF = A1FREQ,Effect = BETA,Direction = sign(Effect)) 
  aou_data_eur = make_all_minor_alleles(aou_data_eur,direction_index = c(1))
  
  # meta dataset
  ad_by_proxy_data = vroom("working/bell_adbyproxy_intersect.txt",show_col_types = F,col_names = F) %>%
    rename(ID = X1,P = X10)
  hwe_stats = vroom("bell_anc_hwe_midp.txt",show_col_types = F) %>%
    mutate(HWE_All_Flag = MIDP_ALL <= 1e-15,
           HWE_Single_Flag = MIDP_AFR <= 1e-15 | MIDP_AMR <= 1e-15 | MIDP_EAS <= 1e-15 |
             MIDP_EUR <= 1e-15 | MIDP_MID <= 1e-15 | MIDP_SAS <= 1e-15)
  
  for (row in 1:nrow(out_df)) {
    curr_row = out_df[row,]
    if (curr_row$ID %in% ad_data$ID | curr_row$IDrev %in% ad_data$ID) {
      match = which(ad_data$ID == curr_row$ID | ad_data$ID == curr_row$IDrev)
      if (length(match)>1) match = which(ad_data$ID == curr_row$ID)
      out_df$AD_P[[row]] = ad_data[[match,'P']]
    }
    if (curr_row$ID %in% ukb_data$ID | curr_row$IDrev %in% ukb_data$ID)
      out_df$UKB_P[[row]] = ukb_data[[which(ukb_data$ID == curr_row$ID | ukb_data$ID == curr_row$IDrev),'P']]
    # aou all participant
    if (curr_row$ID %in% aou_data$ID | curr_row$IDrev %in% aou_data$ID) {
      out_df$AOU_OR[[row]] = exp(aou_data[[which(aou_data$ID == curr_row$ID | aou_data$ID == curr_row$IDrev),'Effect']])
      out_df$AOU_P[[row]] = aou_data[[which(aou_data$ID == curr_row$ID | aou_data$ID == curr_row$IDrev),'P']]
    }
    
    # aou amr
    if (curr_row$ID %in% aou_data_amr$ID | curr_row$IDrev %in% aou_data_amr$ID) {
      out_df$AOU_AMR_P[[row]] = aou_data_amr[[which(aou_data_amr$ID == curr_row$ID | aou_data_amr$ID == curr_row$IDrev),'P']]
    }
    
    # aou afr
    if (curr_row$ID %in% aou_data_afr$ID | curr_row$IDrev %in% aou_data_afr$ID) {
      out_df$AOU_AFR_P[[row]] = aou_data_afr[[which(aou_data_afr$ID == curr_row$ID | aou_data_afr$ID == curr_row$IDrev),'P']]
    }
    
    # aou eur
    if (curr_row$ID %in% aou_data_eur$ID | curr_row$IDrev %in% aou_data_eur$ID) {
      out_df$AOU_EUR_P[[row]] = aou_data_eur[[which(aou_data_eur$ID == curr_row$ID | aou_data_eur$ID == curr_row$IDrev),'P']]
    }
    if (curr_row$ID %in% ad_by_proxy_data$ID | curr_row$IDrev %in% ad_by_proxy_data$ID)
      out_df$AD_by_proxy_P[[row]] = ad_by_proxy_data[[which(ad_by_proxy_data$ID == curr_row$ID | ad_by_proxy_data$ID == curr_row$IDrev),'P']]
    
    # for HWE columns
    if (curr_row$ID %in% hwe_stats$ID | curr_row$IDrev %in% hwe_stats$ID) {
      out_df$HWE_ALL[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_ALL']]
      out_df$HWE_AFR[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_AFR']]
      out_df$HWE_AMR[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_AMR']]
      out_df$HWE_EAS[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_EAS']]
      out_df$HWE_EUR[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_EUR']]
      out_df$HWE_MID[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_MID']]
      out_df$HWE_SAS[[row]] = hwe_stats[[which(hwe_stats$ID == curr_row$ID | hwe_stats$ID == curr_row$IDrev),'MIDP_SAS']]
    }
  }
  
  # Add power values
  power = write_power_df() %>% select(ID,Power_AOU)
  out_df = merge(out_df,power,by='ID') %>% relocate(Power_AOU,.after=AOU_P)
  
  return(out_df %>% filter(!is.na(Bell_P)))
}
getClinicalRiskScore = function() {
  # use Bellenguez hits
  # Note that a number of Bellenguez 
  grs_hits = read_xlsx("Bellenguez_2022_SuppTables.xlsx",sheet = "Supplementary Table 32",skip=2) %>%
    arrange(Chromosome,Position)
  grs_df = data.frame(ID=glue("{grs_hits$Chromosome}:{grs_hits$Position}:{grs_hits$`Other allele`},{grs_hits$`Risk allele`}"),
                      EA=grs_hits$`Risk allele`,Effect=grs_hits$Beta,Pval=grs_hits$`P value`) 
  plink_pos_filter = data.frame(CHROM=grs_hits$Chromosome,POSSTART=grs_hits$Position,POSSTOP=grs_hits$Position)
  
  vroom_write(plink_pos_filter,"grs_hits_plink_pos.bed")
  vroom_write(grs_df,"grs_df_for_score.txt")
}
make_files_for_hwe = function() {
  # get CHRPOS of key variants for HWE testing. Not using IDS here as it is slower
  # The output also includes IDs, so I can do matching there.
  nia_hits = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_p_1e-5.txt",show_col_types = F)
  meta_nia_nimh_hits = vroom("final_meta_results/meta_niagads_nimh_nohwe_chrposrefalt_p_1e-5.TBL",show_col_types = F) %>%
    mutate(CHR = ifelse(CHR == "X","23",CHR),CHR = as.numeric(CHR))
  meta_nia_nimh_unmatched_hits = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/final_meta_results/meta_nia_nimh_nonmatching_p_1e-5_chrposrefalt.TBL")
  ukb_hits = vroom("../Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_p_1e-5.regenie",show_col_types = F)
  aou_hits = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_p_1e-5.txt",show_col_types = F)
  meta_ukb_aou = vroom("final_meta_results/meta_ukb_aou_chrposrefalt_p_1e-5.TBL",show_col_types = F)
  meta_ukb_aou_unmatched_hits = vroom("final_meta_results/meta_ukb_aou_nonmatching_chrposrefalt_p_1e-5.TBL")
  
  bed_file = data.frame(CHR = nia_hits$chr,POSS=nia_hits$pos,POSE=nia_hits$pos) %>%
    add_row(CHR=meta_nia_nimh_hits$CHR,POSS=meta_nia_nimh_hits$POS,POSE=meta_nia_nimh_hits$POS) %>%
    add_row(CHR=meta_nia_nimh_unmatched_hits$CHR,POSS=meta_nia_nimh_unmatched_hits$POS,POSE=meta_nia_nimh_unmatched_hits$POS) %>%
    add_row(CHR=ukb_hits$`#CHROM`,POSS=ukb_hits$GENPOS,POSE=ukb_hits$GENPOS) %>%
    add_row(CHR=aou_hits$CHROM,POSS=aou_hits$GENPOS,POSE=aou_hits$GENPOS) %>%
    add_row(CHR=meta_ukb_aou$CHR,POSS=meta_ukb_aou$POS,POSE=meta_ukb_aou$POS) %>%
    add_row(CHR=meta_ukb_aou_unmatched_hits$CHR,POSS=meta_ukb_aou_unmatched_hits$POS,POSE=meta_ukb_aou_unmatched_hits$POS) %>%
    distinct(.keep_all = T)

  vroom_write(bed_file %>% filter(CHR != 23),"working/all_study_hits.bed")
  print("Wrote BED file for AoU HWE testing to: working/all_study_hits.bed")
  print("Refer to this for code: https://github.com/juliandwillett/AoU_AlzheimersGWAS/blob/main/1.AoU_GWAS/K.CheckForFalsePositives_Anc.sh")
  # print("When running on chr23 in AoU, remember to recode the SEX column in the \\
  #       PSAM file to have the actual sex, given how number of copies can vary with sex")
  return(bed_file)
}
PermDirectionQC = function(df) {
  out_df = df 
  hwe_aou_data = vroom("anc_hwe_midp.txt",show_col_types = F)
  
  rows_to_drop = numeric()
  for (row in 1:nrow(out_df)) {
    r = out_df[row,]
    
    # test by seeing if HWE p val is significant in a single-ancestry cohort, ie deviating
    hwe_match = hwe_aou_data %>% filter(ID == r$ID | IDrev == r$ID) %>%
      filter(MIDP_EUR <= 1e-15 | MIDP_AFR <= 1e-15 | MIDP_AMR <= 1e-15 |
               MIDP_SAS <= 1e-15 | MIDP_EAS <= 1e-15 | MIDP_MID <= 1e-15)
    if (nrow(hwe_match) > 0) rows_to_drop %<>% append(row) 
    else if (nrow(hwe_aou_data %>% filter(ID == r$ID | IDrev == r$ID)) < 1)
      print(glue("Missing HWE for AoU: {r$IDfor}"))
  }
  out_df = out_df[-rows_to_drop,]
  out_df_same_dir = out_df %>% filter(!(str_detect(Direction,"\\+") & str_detect(Direction,"\\-")))
  print(glue("Rows dropped due to HWE: {length(rows_to_drop)}"))
  print(glue("Rows with varying direction of effect: {nrow(out_df) - nrow(out_df_same_dir)}"))
  
  return(out_df)
}
getNumberCommonRareLoci = function(df,maf_cut = 0.01,return_lead,cohort) {
  num_common = num_rare = 0
  num_new_common = num_new_rare = 0
  
  out_df = df %>% filter(Locus == -1) # so get empty data frame to add to for returning lead
  for (loc in sort(unique(df$Locus))) {
    if (loc %% 50 == 0) print(glue("On locus {loc} of {length(unique(df$Locus))}"))
    tmp = df %>% filter(Locus == loc) # all are GW significant
    
    if (return_lead) {
      lead = tmp %>% filter(!!sym(glue("{cohort}_P")) == min(tmp[[glue("{cohort}_P")]]))
      out_df %<>% add_row(lead)
      if (max(lead[[glue("{cohort}_MAF")]]) >= maf_cut & max(lead[[glue("{cohort}_MAF")]]) < (1-maf_cut)) { # testing lead variant
        num_common = num_common + 1
        if (lead$NewOld[[1]] == "New") num_new_common = num_new_common + 1
      }else {
        num_rare = num_rare + 1
        if (lead$NewOld[[1]] == "New") num_new_rare = num_new_rare + 1
      }
    }else{
      if (nrow(tmp %>% filter(!!sym(glue("{cohort}_MAF")) >= maf_cut)) > 0) { # testing all variants in locus
        num_common = num_common + 1
        if (TRUE %in% str_detect(tmp$NewOld,"New")) num_new_common = num_new_common + 1
      }else {
        num_rare = num_rare + 1
        if (TRUE %in% str_detect(tmp$NewOld,"New")) num_new_rare = num_new_rare + 1
      }
    }
  }
  out_df = get_rsids_only(out_df)
  
  print(glue("Total num variants: {nrow(df)}"))
  print(glue("Num common loci: {num_common}. Num rare loci: {num_rare}"))
  print(glue("Num new common loci: {num_new_common}. Num new rare loci: {num_new_rare}"))
  print("Common locus: has at least one common variant in the locus, which could not be the lead")
  print("Rare locus: has ALL rare variants in the locus, including the lead")
  return(out_df)
}
make_all_minor_alleles = function(df,direction_index) {
  out_df = df %>% mutate(Flipped = F)
  for (row in 1:nrow(out_df)) {
    if (out_df$MAF[[row]] >= 0.5) {
      out_df$MAF[[row]] = 1 - out_df$MAF[[row]]
      out_df$Flipped[[row]] = T
      out_df$Effect[[row]] = out_df$Effect[[row]] * -1
      split_id = str_split(out_df$ID[[row]],pattern = "-")[[1]]
      out_df$ID[[row]] = glue("{split_id[[1]]}-{split_id[[2]]}-{split_id[[4]]}-{split_id[[3]]}")
      for (index in direction_index) {
        char = str_sub(out_df$Direction[[row]],index,index)
        if (str_detect(char,"\\+"))
          str_sub(out_df$Direction[[row]],index,index) = "-"
        else if (str_detect(char,"\\-"))
          str_sub(out_df$Direction[[row]],index,index) = "+"
      }
    }
  }
  return(out_df)
}
make_table_lead_loci = function(dataset,in_df,maf_cut = 0.01) {
  out_df = data.frame()

  # Files for getting Z scores and checking results
  # nia_int = vroom("working/niagads_intersects_chrpos.txt",show_col_types = F)
  # nimh_int = vroom("working/nimh_intersects_chrpos.txt",show_col_types = F)
  # nia_nimh_meta_int = vroom("working/niagads_nimh_meta_intersects_chrpos.txt",show_col_types = F)
  # ukb_int = vroom("working/ukb_intersects_chrpos.txt",show_col_types = F)
  # aou_int = vroom("working/aou_intersects_chrpos.txt",show_col_types = F)
  # ukb_aou_meta_int = vroom("working/ukb_aou_meta_intersects_chrpos.txt",show_col_types = F)
  
  if (dataset == "NIA_NIMH") {
    out_df = in_df %>% arrange(CHR,POS) %>% mutate(Effect = NA) %>%
      # First do NIA_NIMH
      mutate(Direction = NIA_NIMH_DIRECTION,MAF = NIA_NIMH_META_MAF) %>%
      make_all_minor_alleles(direction_index=c(1:2)) %>%
      mutate(NIA_NIMH_DIRECTION = Direction,NIA_NIMH_META_MAF = MAF) %>%
      # Second do UKB_AOU
      mutate(Direction = UKB_AOU_DIRECTION,MAF = UKB_AOU_META_MAF) %>%
      make_all_minor_alleles(direction_index=c(1:2)) %>%
      mutate(UKB_AOU_DIRECTION = Direction,UKB_AOU_META_MAF = MAF) %>%
      # Then make the master columns the ones that match here
      mutate(MAF = NIA_NIMH_META_MAF) %>% 
      mutate(Direction = paste0(NIA_NIMH_DIRECTION,UKB_AOU_DIRECTION)) %>%
      select(ID,CHR,Rsid,Gene,MAF,EA,Direction,NIA_P,NIMH_P,NIA_NIMH_META_P,UKB_AOU_META_P,NewOld,IDrev,UKB_P,AOU_P,POS,Flipped) %>%
      mutate(UKB_AOU_Z = NA,.before = UKB_AOU_META_P)
    out_df = get_rsids_only(out_df)
    
    # Get the UKB-AOU Z SCORES
    ukb_aou_meta_int = vroom("working/ukb_aou_meta_intersects_chrpos.txt",show_col_types = F) %>%
      rename(MAF = Freq1,ID=MarkerName,EA=Allele1,NEA=Allele2) %>%
      make_all_minor_alleles(direction_index = c(1:2)) 
    for (row in 1:nrow(out_df)) {
      print(row)
      
      proxy_match = ukb_aou_meta_int %>% filter(ID == out_df$ID[[row]] | ID == out_df$IDrev[[row]])
      if (nrow(ukb_aou_meta_int %>% filter(ID == out_df$ID[[row]],EA != out_df$EA[[row]])) == 0 & nrow(proxy_match) > 0)
        print(paste(out_df$ID[[row]],'varying effect allele'))
      
      if (nrow(proxy_match) > 0) {
        out_df$UKB_AOU_Z[[row]] = round(proxy_match$Effect / proxy_match$StdErr,3)
        out_df$Direction[[row]] = paste0(str_sub(out_df$Direction[[row]],1,2),proxy_match$Direction)
      }
    }
    
    # Manual check. Somehow was not marked in meta-analysis when present, so perhaps single error
    out_df %<>% mutate(UKB_AOU_Z = ifelse(ID == "19-44895007-C-T",-7.37,UKB_AOU_Z)) %>%
      mutate(UKB_AOU_Z = ifelse(ID == "19-44864520-C-T",-1.75,UKB_AOU_Z)) %>%
      mutate(NIMH_P = ifelse(Rsid == "rs157588",0.00929,NIMH_P)) %>%
      mutate(Direction = ifelse(ID == "19-44895007-C-T","----",Direction)) %>%
      mutate(NIA_NIMH_META_P = ifelse(ID == "19-44895007-C-T",6.81e-18,NIA_NIMH_META_P))
  }else if (dataset == "UKB_AOU") {
    # Issue is that some alleles are still major in NIA_NIMH meta and do not match
    # So flip each set sequentially, then make direction column a merged version
    
    out_df = in_df %>% arrange(CHR,POS) %>% mutate(Effect = NA) %>%
      # First do NIA_NIMH
      mutate(Direction = NIA_NIMH_DIRECTION,MAF = NIA_NIMH_META_MAF) %>%
      make_all_minor_alleles(direction_index=c(1:2)) %>%
      mutate(NIA_NIMH_DIRECTION = Direction,NIA_NIMH_META_MAF = MAF) %>%
      # Second do UKB_AOU
      mutate(Direction = UKB_AOU_DIRECTION,MAF = UKB_AOU_META_MAF) %>%
      make_all_minor_alleles(direction_index=c(1:2)) %>%
      mutate(UKB_AOU_DIRECTION = Direction,UKB_AOU_META_MAF = MAF) %>%
      # Then make the master columns the ones that match here
      mutate(MAF = UKB_AOU_META_MAF) %>% 
      mutate(Direction = paste0(NIA_NIMH_DIRECTION,UKB_AOU_DIRECTION)) %>%
      select(ID,CHR,Rsid,Gene,MAF,EA,Direction,NIA_NIMH_META_P,UKB_P,AOU_P,UKB_AOU_META_P,NewOld,AOU_LESS20,IDrev,POS,Flipped) %>%
      mutate(Clin_AD_Z = NA,.before = NIA_NIMH_META_P) %>%
      mutate(EA = ifelse(ID == "19-44838691-A-G",'A',EA))
    out_df = get_rsids_only(out_df) %>% mutate(MAF = round(MAF,3))

    # Get Z score
    nia_nimh_meta_int = vroom("working/niagads_nimh_meta_intersects_chrpos.txt",show_col_types = F) %>%
      rename(MAF = Freq1,ID=MarkerName,EA=Allele1,NEA=Allele2) %>% 
      mutate(Effect = NA) %>%
      make_all_minor_alleles(direction_index = c(1:2)) %>%
      mutate(Zscore = ifelse(Flipped,Zscore * -1,Zscore))
    for (row in 1:nrow(out_df)) {
      print(row)
      
      clin_match = nia_nimh_meta_int %>% filter(ID == out_df$ID[[row]] | ID == out_df$IDrev[[row]])
      if (nrow(nia_nimh_meta_int %>% filter(ID == out_df$ID[[row]],EA != out_df$EA[[row]])) == 0 & nrow(clin_match) > 0)
        print(paste(out_df$ID[[row]],'varying effect allele'))
      
      if (nrow(clin_match) > 0) out_df$Clin_AD_Z[[row]] = round(clin_match$Zscore,3)
    }
    return(out_df)
  }

  return(out_df %>% 
           mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
                  Gene = str_replace_all(Gene," ",""),
                  Gene = gsub("\\([^)]*\\)", "", Gene),
                  Gene = str_replace_all(Gene,",",", ")))
}
makeTable3_4_CADD_Hits_with_FAVOR = function(df,cadd_thresh=10,maf_cut = 0.01) {
# Organize results
  favor_annot = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",col_names=F,delim=",",skip=1)
  names(favor_annot) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",n_max=1,delim="\t"))
  df_phred = df %>% #filter(ID %in% (favor_annot %>% filter(cadd_phred >= cadd_thresh))$variant_vcf) %>%
    filter(NewOld == "New",EnhancerRole == "Yes" | EpiActivMax >= cadd_thresh | EpiReprMax >= cadd_thresh |
             EpiTransMax >= cadd_thresh | CADD_Phred >= cadd_thresh) 
  
  # Remove rare variants that are not in both datasets (for immediate QC on single-variant GWAS)
  df_arranged = df_phred %>% filter(A1Freq >= maf_cut) %>% add_row(df_phred %>% filter(A1Freq < maf_cut))

  table3 = df_arranged %>% select(ID,CHROM,Rsid,ProximalGene,A1Freq,Direction,OR,CI,MetaP,AoU_Pval,UKB_Pval,NewOld) 
  table4 = df_arranged %>% select(ID,CHROM,Rsid,ProximalGene,CADD_Phred,EnhancerRole,EpiActivMax,EpiReprMax,
                               EpiTransMax,Conserved)
  
  write_xlsx(list("Table3.CADD_Variant_GWAS_Info" = table3,"Table4.CADD_Variant_FAVOR_Scores" = table4),
             path = "Paper_Tables_Figures/Table3_4_CADD_Variants_FAVOR.xlsx")
  print("Wrote table to: Paper_Tables_Figures/Paper_Tables_Figures/Table3_4_CADD_Variants_FAVOR.xlsx")
  
  return(list(table3,table4))
}
add_cell_expression = function(df,p_cutoff = 0.01) {
  print("Focus on clarifying temporal variables, ie comparing mild AD outcome to null outcome (second to first), then severe to mild (third to second)")
  print("These data should represent data from samples cleaned with DecontX")
  
  # out_df = data.frame(Outcome=character(),Comparison=character(),Gene=character(),
  #                     CellPop=character(),AvgLog2FC=numeric(),PValAdj=numeric())
  # for (f in list.files("~/true_lab_storage/00_AoU/mathys_deg_data_precomputed/ROSMAP_snRNAseq_PFC/Results/Disease_progression",recursive=T,pattern=".csv",full.names = T)) {
  #   spl_str = str_split(f,"/")[[1]]
  #   cell_data = read_csv(f)
  #   out_df %<>% add_row(Outcome = spl_str[[11]],Comparison=spl_str[[13]],
  #                       Gene = cell_data$gene,CellPop=cell_data$cluster_id,
  #                       AvgLog2FC = cell_data$logFC,PValAdj=cell_data$p_adj.glb)
  # }
  # vroom_write(out_df,glue("scrna_logfold_stats_by_cell.txt"))
  
  sc_results = vroom("scrna_logfold_stats_by_cell.txt") %>%
    filter(PValAdj <= p_cutoff) %>% arrange(Gene,desc(AvgLog2FC)) %>% 
    mutate(Comparison = str_replace(Comparison,".csv","")) %>% group_by(Gene) %>%
    mutate(CellPop = ifelse(str_detect(Comparison,"second_vs_first"),glue("{CellPop} 2v1"),CellPop)) %>%
    mutate(CellPop = ifelse(str_detect(Comparison,"third_vs_second"),glue("{CellPop} 3v2"),CellPop))

  out_df = df %>% mutate(AmyloidProgr_Cells_P_LogFC = NA,GpathProgr_Cells_P_LogFC=NA,
                         NFTProgr_Cells_P_LogFC=NA,PlaqDProgr_Cells_P_LogFC=NA,
                         PlaqNProgr_Cells_P_LogFC=NA,TanglesProgr_Cells_P_LogFC=NA)
  for (row in 1:nrow(out_df)) {
    if (row %% 100 == 0) print(glue("On row {row} of {nrow(out_df)}"))
    row_genes = out_df[row,] %>% select(Gene,EnhancerLinkedGene) %>%
      mutate(Gene = paste0(Gene,",",EnhancerLinkedGene,collapse=",")) %>%
      mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
             Gene = str_replace_all(Gene," ",""),
             Gene = gsub("\\([^)]*\\)", "", Gene)) %>%
      separate_rows(Gene, sep = ",") %>% distinct()
    
    for (g in unique(row_genes$Gene)) { # cross reference gene, unique to avoid duplicates
      sc_res_match = sc_results %>% filter(Gene == g)
      amyloid_match = sc_res_match %>% filter(Outcome == "amyloid") 
      gpath_match = sc_res_match %>% filter(Outcome == "gpath") 
      nft_match = sc_res_match %>% filter(Outcome == "nft") 
      plaqd_match = sc_res_match %>% filter(Outcome == "plaq_d")
      plaqn_match = sc_res_match %>% filter(Outcome == "plaq_n")
      tangles_match = sc_res_match %>% filter(Outcome == "tangles")

      if (nrow(amyloid_match)>0) out_df$AmyloidProgr_Cells_P_LogFC[[row]] = 
        paste(out_df$AmyloidProgr_Cells_P_LogFC[[row]],paste(glue("({g})"),
                    glue("({amyloid_match$CellPop};{amyloid_match$PValAdj};{amyloid_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(gpath_match)>0) out_df$GpathProgr_Cells_P_LogFC[[row]] = 
        paste(out_df$GpathProgr_Cells_P_LogFC[[row]],paste(glue("({g})"),
                    glue("({gpath_match$CellPop};{gpath_match$PValAdj};{gpath_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(nft_match)>0) out_df$NFTProgr_Cells_P_LogFC[[row]] = 
        paste(out_df$NFTProgr_Cells_P_LogFC[[row]],paste(glue("({g})"),
                    glue("({nft_match$CellPop};{nft_match$PValAdj};{nft_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(plaqd_match)>0) out_df$PlaqDProgr_Cells_P_LogFC[[row]] = 
        paste(out_df$PlaqDProgr_Cells_P_LogFC[[row]],paste(glue("({g})"),
                    glue("({plaqd_match$CellPop};{plaqd_match$PValAdj};{plaqd_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(plaqn_match)>0) out_df$PlaqNProgr_Cells_P_LogFC[[row]] = 
        paste(out_df$PlaqNProgr_Cells_P_LogFC[[row]],paste(glue("({g})"),
                    glue("({plaqn_match$CellPop};{plaqn_match$PValAdj};{plaqn_match$AvgLog2FC})"),collapse=","),",")
      if (nrow(tangles_match)>0) out_df$TanglesProgr_Cells_P_LogFC[[row]] = 
        paste(out_df$TanglesProgr_Cells_P_LogFC[[row]],paste(glue("({g})"),
                    glue("({tangles_match$CellPop};{tangles_match$PValAdj};{tangles_match$AvgLog2FC})"),collapse=","),",")
    }
  }
  # clean up cell name columns a bit
  for (col in 29:34) {
    out_df[[col]] = str_replace(out_df[[col]],"NA ","")
    out_df[[col]] = str_replace_all(out_df[[col]]," ","_")
  }
  return(out_df)
}
get_logfc_sign = function(v,ce,co) { # take in vector of strings, cell, and comparison
  lines = v[which(str_detect(v,ce) & str_detect(v,co))]
  if (length(lines) < 1) return("")
  out_str = ""
  for (line in lines) {
    spl_l = str_split(line,";")[[1]][[3]]
    if (str_detect(spl_l,"-")) out_str = paste0(out_str,"-")
    else out_str = paste0(out_str,"+")
  }
  return(out_str)
}
make_table_cell_expression = function(df,col,clean_cell_columns,include_signs,only_new,enhancer_role) {
  filt_df = df %>% 
    mutate(OutcomeCellSignificant = str_replace(OutcomeCellSignificant,"NA","")) %>% 
    drop_na(OutcomeCellSignificant)
  if (enhancer_role) filt_df %<>% filter(EnhancerRole == "Yes")
  
  # Make df
  out_df = data.frame(ID=filt_df$ID,CHR=filt_df$CHR,Rsid=filt_df$Rsid,MAF=filt_df[[paste0(col,"_MAF")]],
                      Direction = paste0(filt_df$NIA_NIMH_DIRECTION,filt_df$UKB_AOU_DIRECTION),
                      Pval=filt_df[[paste0(col,"_P")]],Locus=filt_df$Locus,
                      Gene=filt_df$Gene,LinkedGene=filt_df$EnhancerLinkedGene,
                      Cell=filt_df$OutcomeCellSignificant,Locus=filt_df$Locus,NewOld=filt_df$NewOld,
                      UKB_P = filt_df$UKB_P,AOU_P = filt_df$AOU_P,
                      NIA_NIMH_DIRECTION = filt_df$NIA_NIMH_DIRECTION,
                      UKB_AOU_DIRECTION = filt_df$UKB_AOU_DIRECTION) %>%
    mutate(Gene = str_replace_all(Gene,"\\(dist.*?\\)", ""),
           Gene = str_replace_all(Gene," ",""),
           Gene = gsub("\\([^)]*\\)", "", Gene)) %>%
    mutate(Direction = paste0(NIA_NIMH_DIRECTION,UKB_AOU_DIRECTION))
  
  # Flip alleles and plug in rsids
  index = ifelse(col == "NIA_NIMH_META",c(1:2),c(3:4))
  out_df = make_all_minor_alleles(out_df %>% mutate(Effect = NA),index)
  out_df = get_rsids_only(out_df) %>% mutate(MAF = round(MAF,4))
  out_df = match_signs_proxy(out_df)
  
  # Stick in DEG information
  diff_expr_genes = character()
  for (row in 1:nrow(out_df)) {
    matching_locus = (df %>% filter(Locus == out_df$Locus[[row]]))
    genes = unique(str_split(paste(out_df$Gene[[row]],out_df$LinkedGene[[row]],sep=","),",")[[1]])
    
    # Count cell types and give number of subpopulations for each type (to make more concise)
    clean_str = "" # for each row
    if (clean_cell_columns) {
      spl = str_split(out_df$Cell[[row]],",")[[1]]

      for (gene in genes) { # to deal with more than one gene
        # Do work
        gene_edit = str_replace_all(gene,"\\(dist.*?\\)", "")
        gene_edit = str_replace_all(gene_edit," ","")
        gene_edit = gsub("\\([^)]*\\)", "", gene_edit)
        
        spl_g = grep(glue("\\({gene_edit}\\)"),spl,value=T) # focus on entries for given gene

        cells = c("Exc ","Inh ","Oli ","OPC ","Ast ","Mic ")
        comps = c(" 2v1"," 3v1"," 4v1"," 3v2"," 4v2"," 4v3")
        
        for (cell in cells) {
          for (comp in comps) {
            count = length(which(str_detect(spl_g,cell) & str_detect(spl_g,comp)))
            # Also want to pull by p-value and the log2fc (2nd and 3rd entry delim by semicolons)
            if (count > 0) {
              # Sign sometimes vary by population, so less informative, thus include only when selected
              if (include_signs) {
                sign = paste0(get_logfc_sign(spl_g,cell,comp))
              }else sign = ""
              clean_str = glue("{clean_str} \\({gene_edit}\\) {cell} {comp} (x{count}) {sign}")
              
              if (only_new) {
                if (out_df$NewOld[[row]] == "New")
                  diff_expr_genes %<>% append(gene)
              }else
                diff_expr_genes %<>% append(gene)
            }
          }
        }
      }
      
    }else clean_str = out_df[[tabl]]$Cell[[row]] #
    out_df$Cell[[row]] = str_squish(clean_str)
    out_df$Cell[[row]] = gsub("\\\\","",out_df$Cell[[row]])

    # remove gene if there is only one being evaluated
    if (length(genes) == 1) {
      out_df$Cell[[row]] = str_replace_all(out_df$Cell[[row]],glue("\\({genes[[1]]}\\) "),"")
    }
  }
  out_df %<>% filter(Cell != "")
  if (only_new) out_df %<>% filter(NewOld == "New")
  
  # Remove Genes from LinkedGene not in Mathys et al data, to clean it up more
  if (enhancer_role) {
    for (row in 1:nrow(out_df)) {
      curr_row = out_df[row,]
      curr_genes = str_split(curr_row$LinkedGene,pattern = ",")[[1]]
      for (gene in curr_genes) {
        if (!str_detect(curr_row$Cell,gene))
          out_df$LinkedGene[[row]] = str_replace(out_df$LinkedGene[[row]],gene,'')
      }
    }
    out_df %<>% mutate(LinkedGene = str_replace_all(LinkedGene,",,",","))
  }
  # Count genes
  genes = length(unique(diff_expr_genes))

  # Select most significant variant from a locus if multiple variants present
  out_df %<>% group_by(Locus) %>% filter(Pval == min(Pval) | !is.na(LinkedGene)) %>% ungroup(Locus)

  # Report statistics and return count
  print(paste("Number of genes differentially expressed:",genes))
  return(out_df)
}
make_table_cell_expression_bypath = function(df,coh_col,clean_cell_columns,include_signs=T,only_new,
                                                   enhancer_linked) {
  filt_df = df %>% mutate(MAF = !!sym(paste0(coh_col,"_MAF")),P = !!sym(paste0(coh_col,"_P"))) %>%
    mutate(Direction = paste0(NIA_NIMH_DIRECTION,UKB_AOU_DIRECTION))
  if (only_new) filt_df %<>% filter(NewOld == "New")
  
  # Flip alleles and plug in rsids
  index = ifelse(coh_col == "NIA_NIMH_META",c(1:2),c(3:4))
  filt_df = make_all_minor_alleles(filt_df %>% mutate(Effect = NA),direction_index = index)
  filt_df = get_rsids_only(filt_df) %>% mutate(MAF = round(MAF,4)) %>%
    group_by(Locus) %>% filter(P == min(P) | !is.na(EnhancerLinkedGene)) %>% ungroup(Locus)
  
  amy_table = filt_df %>% filter(!is.na(AmyloidProgr_Cells_P_LogFC))
  gpath_table = filt_df %>% filter(!is.na(GpathProgr_Cells_P_LogFC))
  nft_table = filt_df %>% filter(!is.na(NFTProgr_Cells_P_LogFC))
  plaqd_table = filt_df %>% filter(!is.na(PlaqDProgr_Cells_P_LogFC))
  plaqn_table = filt_df %>% filter(!is.na(PlaqNProgr_Cells_P_LogFC))
  tangles_table = filt_df %>% filter(!is.na(TanglesProgr_Cells_P_LogFC))
  tables = list(amy_table,gpath_table,nft_table,plaqd_table,plaqn_table,tangles_table)
  
  out_df = list()
  for (tabl in 1:6) {
    out_df[[tabl]] = data.frame(ID=tables[[tabl]]$ID,CHR=tables[[tabl]]$CHR,
                                Rsid=tables[[tabl]]$Rsid,MAF=tables[[tabl]]$MAF,
                                Gene=tables[[tabl]]$Gene,
                                  LinkedGene=tables[[tabl]]$EnhancerLinkedGene,
                                Cell=tables[[tabl]][[(31+tabl)]],
                                Direction = paste0(tables[[tabl]]$NIA_NIMH_DIRECTION,
                                                   tables[[tabl]]$UKB_AOU_DIRECTION))
    out_df[[tabl]] = out_df[[tabl]] 
    
    for (row in 1:nrow(out_df[[tabl]])) {
      if (enhancer_linked)
        genes = unique(str_split(paste(out_df[[tabl]]$Gene[[row]],out_df[[tabl]]$LinkedGene[[row]],sep=","),",")[[1]])
      else
        genes = unique(str_split(out_df[[tabl]]$Gene[[row]],",")[[1]])
      # Count cell types and give number of subpopulations for each type (to make more concise)
      clean_str = "" # for each row
      if (clean_cell_columns) {
        spl = str_split(out_df[[tabl]]$Cell[[row]],",")[[1]]
        for (gene in genes) { # to deal with more than one gene
          gene_edit = str_replace_all(gene,"\\(dist.*?\\)", "")
          gene_edit = str_replace_all(gene_edit," ","")
          gene_edit = gsub("\\([^)]*\\)", "", gene_edit)

          spl_g = grep(glue("({gene_edit})"),spl,value=T) # focus on entries for given gene
          if (tabl > 3) spl_g = grep(gene_edit,spl,value=T)

          cells = c("Exc","Inh","Oli","OPC","Ast","Mic")
          comps = c("2v1","3v1","3v2")
          
          for (cell in cells) {
            for (comp in comps) {
              count = length(which(str_detect(spl_g,cell) & str_detect(spl_g,comp)))
              if (count > 0) {
                if (include_signs) {
                  sign = paste0(get_logfc_sign(spl_g,cell,comp))
                }else sign = ""
                clean_str = glue("{clean_str} \\({gene_edit}\\) {str_replace(cell,'_','')} {comp} (x{count}) {sign}")
              }
            }
          }
        }
      }else clean_str = out_df[[tabl]]$Cell[[row]] # 
      out_df[[tabl]]$Cell[[row]] = clean_str
      out_df[[tabl]]$Cell[[row]] = gsub("\\\\","",out_df[[tabl]]$Cell[[row]] )
      
      # remove gene name if there is only one being evaluated
      if (length(genes) == 1) {
        out_df[[tabl]]$Cell[[row]] = str_replace_all(out_df[[tabl]]$Cell[[row]],glue("\\({genes[[1]]}\\) "),"")
      }
    }
  }
  names(out_df) = c("Amyloid","GPath","NFT","PlaqD","PlaqN","Tangles")
  out_full_df = out_df[[1]] %>% mutate(Path=names(out_df)[[1]],.before=CHR) %>%
    add_row(out_df[[2]] %>% mutate(Path=names(out_df)[[2]],.before=CHR)) %>%
    add_row(out_df[[3]] %>% mutate(Path=names(out_df)[[3]],.before=CHR)) %>%
    add_row(out_df[[4]] %>% mutate(Path=names(out_df)[[4]],.before=CHR)) %>%
    add_row(out_df[[5]] %>% mutate(Path=names(out_df)[[5]],.before=CHR)) %>%
    add_row(out_df[[6]] %>% mutate(Path=names(out_df)[[6]],.before=CHR))  %>%
    relocate(Path,.after=LinkedGene)
  
  # DO QC
  out_full_df %<>% filter(Cell != "") %>% arrange(CHR,Gene) %>%
    mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
           Gene = str_replace_all(Gene," ",""),
           Gene = gsub("\\([^)]*\\)", "", Gene))
  if (enhancer_linked) out_full_df %<>% filter(!is.na(LinkedGene))
  
  # CLEAN UP LINKED GENE COLUMN
  for (row in 1:nrow(out_full_df)) {
    curr_row = out_full_df[row,]
    if (is.na(curr_row$LinkedGene)) next
    curr_genes = str_split(curr_row$LinkedGene,pattern = ",")[[1]]
    for (gene in curr_genes) {
      if (!str_detect(curr_row$Cell,gene))
        out_full_df$LinkedGene[[row]] = str_replace(out_full_df$LinkedGene[[row]],gene,'')
    }
  }
  out_full_df %<>% mutate(LinkedGene = str_replace_all(LinkedGene,",,",","))
  
  return(out_full_df)
}
# getNumberPhredAnnot = function(df,maf_cut = 0.01) {
#   favor_annot = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",col_names=F,delim=",",skip=1,show_col_types = F)
#   names(favor_annot) = names(vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/favor_hits_annot.txt",n_max=1,delim="\t",show_col_types = F))
#   
#   cadd_hits = df %>% filter(CADD_Phred >= 10)
#   cadd_hits_20 = df %>% filter(CADD_Phred >= 20)
#   num_cadd = c(nrow(cadd_hits %>% filter(A1Freq >= maf_cut)),nrow(cadd_hits %>% filter(A1Freq < maf_cut)))
#   num_cadd_new = c(nrow(cadd_hits %>% filter(A1Freq >= maf_cut,NewOld == "New")),
#                    nrow(cadd_hits %>% filter(A1Freq < maf_cut,NewOld == "New")))
#   print(glue("Num CADD Phred >= 10 Common {num_cadd[[1]]}, Rare {num_cadd[[2]]}"))
#   print(glue("Num New CADD Phred >= 10 Common {num_cadd_new[[1]]}, Rare {num_cadd_new[[2]]}"))
#   print(glue("Num CADD Phred >= 10 and AoU nominal: {nrow(cadd_hits %>% filter(AoU_Pval <= 0.05))}"))
#   num_cadd_20 = c(nrow(cadd_hits_20 %>% filter(A1Freq >= maf_cut)),nrow(cadd_hits_20 %>% filter(A1Freq < maf_cut)))
#   num_cadd_20_new = c(nrow(cadd_hits_20 %>% filter(A1Freq >= maf_cut,NewOld == "New")),
#                       nrow(cadd_hits_20 %>% filter(A1Freq < maf_cut,NewOld == "New")))
#   print(glue("Num CADD Phred >= 20 Common {num_cadd_20[[1]]}, Rare {num_cadd_20[[2]]}"))
#   print(glue("Num New CADD Phred >= 20 Common {num_cadd_20_new[[1]]}, Rare {num_cadd_20_new[[2]]}"))
#   
#   # Print out the rare novel variants
#   print(glue("Novel variants: {cadd_hits %>% filter(NewOld == 'New') %>% select(ID,Rsid,ProximalGene,CADD_Phred)}"))
#   
#   ###############
#   # Epi stuff
#   epi_hits = (favor_annot %>% filter(variant_vcf %in% df$ID) %>%
#          filter(encodeh3k4me1_sum > 10 | encodeh3k4me2_sum > 10 | encodeh3k4me3_sum > 10 | 
#                   encodeh3k9ac_sum > 10 | encodeh3k9me3_sum > 10 | encodeh3k27ac_sum > 10 |
#                   encodeh3k27me3_sum > 10 | encodeh3k36me3_sum > 10 | encodeh3k79me2_sum > 10 | 
#                   encodeh4k20me1_sum > 10 | encodeh2afz_sum > 10 | apc_epigenetics_active > 10 | 
#                   apc_epigenetics_repressed > 10 | apc_epigenetics_transcription > 10 | 
#                   apc_transcription_factor > 10))
#   epi_hits_in_df = df %>% filter(ID %in% epi_hits$variant_vcf)
#   num_epi = c(nrow(epi_hits_in_df %>% filter(A1Freq >= maf_cut)),nrow(epi_hits_in_df %>% filter(A1Freq < maf_cut)))
#   print(glue("Num Epi Phred >= 10 Common {num_epi[[1]]}, Rare {num_epi[[2]]}"))
#   print(glue("Num nominal in AoU: {nrow(epi_hits_in_df %>% filter(AoU_Pval <= 0.05))}"))
#   
#   epi_hits_in_df_20 = df %>% filter(ID %in% epi_hits$variant_vcf,EpiActivMax >= 20 |
#                                     EpiReprMax >= 20 | EpiTransMax >= 20)
#   num_epi_20 = c(nrow(epi_hits_in_df_20 %>% filter(A1Freq >= maf_cut)),
#                  nrow(epi_hits_in_df_20 %>% filter(A1Freq < maf_cut)))
#   print(glue("Num Epi Phred >= 20 Common {num_epi_20[[1]]}, Rare {num_epi_20[[2]]}"))
# }
makeFigure5 = function(df,maf_cut = 0.01) {
  # get entries
  all_df = df
  all_count = c(nrow(all_df),nrow(all_df %>% filter(A1Freq>=maf_cut)),nrow(all_df %>% filter(A1Freq<maf_cut)))

  cadd_df = df %>% filter(CADD_Phred >= 10)
  cadd_count = c(nrow(cadd_df),nrow(cadd_df %>% filter(A1Freq>=maf_cut)),nrow(cadd_df %>% filter(A1Freq<maf_cut)))
  
  epi_df = df %>% filter(EpiActivMax >= 10 | EpiReprMax >= 10 | EpiTransMax >= 10)
  epi_count = c(nrow(epi_df),nrow(epi_df %>% filter(A1Freq>=maf_cut)),nrow(epi_df %>% filter(A1Freq<maf_cut)))
  
  # make figure
  fig_df = data.frame(Group=character(),Count=numeric(),Subgroup=character())
  fig_df %<>% add_row(Group="All",Count=all_count[[1]],Subgroup="Total") %>%
    add_row(Group="All",Count=all_count[[2]],Subgroup="Common") %>%
    add_row(Group="All",Count=all_count[[3]],Subgroup="Rare") %>%
    add_row(Group="CADD",Count=cadd_count[[1]],Subgroup="Total") %>%
    add_row(Group="CADD",Count=cadd_count[[2]],Subgroup="Common") %>%
    add_row(Group="CADD",Count=cadd_count[[3]],Subgroup="Rare") %>%
    add_row(Group="Epigenetic",Count=epi_count[[1]],Subgroup="Total") %>%
    add_row(Group="Epigenetic",Count=epi_count[[2]],Subgroup="Common") %>%
    add_row(Group="Epigenetic",Count=epi_count[[3]],Subgroup="Rare")
  fig_df$Subgroup <- factor(fig_df$Subgroup, levels = c("Total", "Common", "Rare"))
  
  plt = ggplot(fig_df,aes(x=Group,y=Count,fill=Subgroup)) + 
    geom_bar(stat="identity",position=position_dodge()) +
    theme_bw() + theme(text = element_text(size=18)) + xlab("Functional Annotation") +
    ylab("Number of Variants")
  print(plt)
  ggsave(plot=plt,filename="Paper_Tables_Figures/Figure5.png",width=6,height=4)
  print("Made file: Paper_Tables_Figures/Figure5.png")
  
  return(fig_df)
}
compareRegenieWithRelated_Vs_WithoutRelated = function() {
  hits_rel = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt") %>%
    mutate(ID = glue("{CHROM}:{GENPOS}:{ALLELE0},{ALLELE1}"))
  hits_nonrel = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NO_RELATED_GWAS/aou_AD_any_anc_all_gwas_bt_nomcc_norel_gwsig.txt") %>%
    filter(LOG10P >= -log10(5e-8)) %>% arrange(CHROM,GENPOS) %>% mutate(Pval=10^(-LOG10P)) %>%
    mutate(ID = glue("{CHROM}:{GENPOS}:{ALLELE0},{ALLELE1}"))
  df = data.frame(ID=unique(hits_rel$ID %>% append(hits_nonrel$ID)),P_rel=NA,P_nonrel=NA) %>% arrange(ID)
  
  rel_lesssig = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse_1e-3.txt") %>%
    mutate(ID = glue("{CHROM}:{GENPOS}:{ALLELE0},{ALLELE1}"))
  nonrel_lesssig = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NO_RELATED_GWAS/aou_AD_any_anc_all_gwas_bt_nomcc_norel_gwsig.txt") %>%
    mutate(Pval=10^(-LOG10P)) %>% mutate(ID = glue("{CHROM}:{GENPOS}:{ALLELE0},{ALLELE1}"))
  missing_ids = character()
  for (row in 1:nrow(df)) {
    print(row)
    print(df$ID[[row]])
    match_rel = rel_lesssig %>% filter(ID == df$ID[[row]]) %>% distinct()
    match_nonrel = nonrel_lesssig %>% filter(ID == df$ID[[row]]) %>% distinct()
    if (nrow(match_rel)>0) df$P_rel[[row]] = match_rel$Pval
    else missing_ids %<>% append(df$ID[[row]])
    if (nrow(match_nonrel)>0) df$P_nonrel[[row]] = match_nonrel$Pval
    else missing_ids %<>% append(df$ID[[row]])
  }
  print(glue("Missing: {unique(missing_ids)}"))
  return(df)
}
isolateSharedHits = function(df) {
  # None of the flagged variants are significant in the meta-analysis
  # ukb_hwe_flagged = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/EUR_hwe_less_1e-15.hardy") %>%
  #   separate(ID,into=c("CHR","POS"),sep = ":") %>%
  #   mutate(CHRPOS = glue("{CHR}-{POS}")) %>% filter(CHRPOS %in% df$CHRPOS) 
  
  tmp = df %>% filter(UKB_Pval < 1 | is.na(UKB_Pval)) %>% # deal with flagged UKB variants that did not get removed
    group_by(Locus) %>% 
    filter(length(which(!is.na(AoU_Pval) & AoU_Pval <= 0.05 &
      (!is.na(Niagads_Pval) & Niagads_Pval <= 0.05) | 
                          (!is.na(UKB_Pval) & UKB_Pval <= 0.05) |
        (!is.na(UKB_Pval) & !is.na(Niagads_Pval))))>0) %>% # check by locus, not by single variant
    ungroup(Locus) 
  
  getNumberCommonRareLoci(tmp,return_lead = F)
  return(tmp)
}
addDiseaseGroup = function(df,p_cutoff=0.01) {
  # Use similar code to add_cell_expression function, because similarly organized data
  
  print("Applied FDR cutoff of 0.05 to raw data. Rerun this early code if I want to be more lenient")
  # out_df = data.frame(Comparison=character(),Gene=character(),
  #                     CellPop=character(),AvgLog2FC=numeric(),PValAdj=numeric())
  # for (f in list.files("~/true_lab_storage/00_AoU/mathys_deg_data_precomputed/ROSMAP_snRNAseq_PFC/Results/Four-group",recursive=T,pattern=".csv",full.names = T)) {
  #   spl_str = str_split(f,"/")[[1]]
  #   compar = str_replace(spl_str[[12]],paste0(spl_str[[11]],"_group_"),"")
  #   compar = str_replace(compar,".csv","")
  #   cell_data = read_csv(f) %>% filter(p_adj.glb <= p_cutoff)
  #   if (nrow(cell_data) > 0)
  #     out_df %<>% add_row(Comparison=compar,
  #                       Gene = cell_data$gene,CellPop=cell_data$cluster_id,
  #                       AvgLog2FC = cell_data$logFC,PValAdj=cell_data$p_adj.glb)
  # }
  # vroom_write(out_df,glue("scrna_logfold_stats_by_cell_group_by_pathology_and_symptoms_1234.txt"))
  
  sc_results = vroom("scrna_logfold_stats_by_cell_group_by_pathology_and_symptoms_1234.txt") %>%
    filter(PValAdj <= p_cutoff) %>% arrange(Gene,desc(AvgLog2FC)) %>% group_by(Gene) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"second_vs_first"),"2v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"third_vs_first"),"3v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_first"),"4v1",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"third_vs_second"),"3v2",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_second"),"4v2",Comparison)) %>%
    mutate(Comparison = ifelse(str_detect(Comparison,"fourth_vs_third"),"4v3",Comparison)) 

  out_df = df %>% mutate(OutcomeCellSignificant = NA)
  for (row in 1:nrow(out_df)) {
    if (row %% 100 == 0) print(glue("On row {row} of {nrow(out_df)}"))
    row_genes = out_df[row,] %>% select(Gene,EnhancerLinkedGene) %>%
      mutate(Gene = paste0(Gene,",",EnhancerLinkedGene,collapse=",")) %>%
      mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
             Gene = str_replace_all(Gene," ",""),
             Gene = gsub("\\([^)]*\\)", "", Gene)) %>%
      separate_rows(Gene, sep = ",")
    
    for (g in unique(row_genes$Gene)) { # cross reference gene, unique to avoid duplicates
      sc_res_match = sc_results %>% filter(Gene == g) %>% mutate(CellPopComp = paste(CellPop,Comparison))
      
      if (nrow(sc_res_match)>0) {
        out_df$OutcomeCellSignificant[[row]] = 
          paste0(out_df$OutcomeCellSignificant[[row]],paste(glue("({g})"),
                glue("({sc_res_match$CellPopComp};{sc_res_match$PValAdj};{sc_res_match$AvgLog2FC})"),
                                                           collapse=","),collapse=",")
        if (length(unique(row_genes$Gene)) > 1 & 
            g != unique(row_genes$Gene)[[length(unique(row_genes$Gene))]])
          out_df$OutcomeCellSignificant[[row]] = paste0(out_df$OutcomeCellSignificant[[row]],",")
      }
    }
  }
  # Clean up text a bit
  out_df$OutcomeCellSignificant = str_replace(out_df$OutcomeCellSignificant,"NA ","")
  
  return(out_df)
}
prep_data_for_separate_loci = function(df,study_str) {
  # Run it by chromosome, so make no assumptions about loci by position.
  out_list = list()
  for (chr in 1:23) {
    out_list[[chr]] = df %>% filter(CHR == chr)
  }
  out_path = glue("locus_hits/locus_hits_{study_str}")
  if (!file.exists(out_path)) dir.create(out_path)
  for (chr in 1:23) {
    vroom_write(out_list[[chr]],glue("{out_path}/chr{chr}_for_clump.txt"))
  }
  print("Files written to: locus_hits")
}
prep_data_for_separate_loci_allcohorts = function(gw) {
  for (coh in c("NIA","NIMH","NIA_NIMH_META","UKB","AOU","UKB_AOU_META")) {
    col_match = glue("{coh}_P")
    df = gw[which(gw[[col_match]] < 5e-8),]
    prep_data_for_separate_loci(df,study_str=coh)
  }
  print("Pipe data into AoU using: https://github.com/juliandwillett/AoU_AlzheimersGWAS/blob/main/1.AoU_GWAS/O.Get_Independent_Loci.sh")
}
separate_loci_gwas = function(df,favor_annot) {
  out_list = list()
  print("Direction differences dealt with here, but not nominal")
  for (coh in c("NIA","NIMH","NIA_NIMH_META","UKB","AOU","UKB_AOU_META")) {
    print(glue("On cohort {coh}"))
    if (coh != "NIMH")
      clumps = vroom(glue("clumps/clumps_{coh}.txt"),show_col_types = F)
    else
      clumps = vroom("clumps/clumps_NIA.txt",show_col_types = F) %>% filter(ID == "")
    curr_df = df[which(df[[glue("{coh}_P")]] <= 5e-8),] %>% 
      mutate(Gene=NA,Rsid=NA,Locus=NA,EnhancerRole=NA,EnhancerLinkedGene=NA,
             ClumpedLocus=FALSE,NewOld=NA,ClinNom=NA)

    # Then start the clumping. Go through curr_df and match to clump row, with that being clump number
    failed_clumping_loci = numeric()
    for (row in 1:nrow(curr_df)) {
      curr_row = curr_df[row,]
      
      # DO CLUMPING
      clump_match = clumps %>% filter(ID == curr_row$ID | ID == curr_row$IDrev |
                                        str_detect(SP2,curr_row$ID) |
                                        str_detect(SP2,curr_row$IDrev)) 
      if (nrow(clump_match) > 1) 
        clump_match = clump_match %>% filter(ID == curr_row$ID | str_detect(SP2,curr_row$ID))
      if (nrow(clump_match) > 0) {
        clump_num = which(clumps$ID == clump_match$ID)
        curr_df$Locus[[row]] = clump_num
        curr_df$ClumpedLocus[[row]] = TRUE
      }
      
      # INCLUDE NEW/OLD DELINEATION
      if (str_detect(coh,"META")) true_coh = str_replace(coh,"_META","")
      else true_coh = coh
      annot_match = favor_annot[[true_coh]] %>% filter(ID == curr_row$ID | ID == curr_row$IDrev) 
      curr_df$NewOld[[row]] = annot_match$NewOld[[1]]
      
      # INCLUDE PROXIMAL GENE
      curr_df$Gene[[row]] = annot_match$ProximalGene[[1]]
      curr_df$EnhancerRole[[row]] = annot_match$EnhancerRole[[1]]
      curr_df$EnhancerLinkedGene[[row]] = annot_match$EnhancerLinkedGene[[1]]
      
      # INCLUDE RSID
      curr_df$Rsid[[row]] = annot_match$Rsid[[1]]
      
      # INCLUDE NUMBER OF DATASETS IT IS NOMINAL IN
      curr_df$ClinNom[[row]] = curr_row$NIA_NIMH_META_P <= 0.05
      
      # INCLUDE IF PASSED HWE
      curr_df$HWE_PASS[[row]] = curr_row$HWE_PASS
    }
    
    # Deal with code artifact in UKB
    if (coh == "UKB") curr_df %<>% filter(UKB_P != 0)
    
    # Then return data, putting clumped variants at the top, removing duplicates
    out_list[[coh]] = curr_df %>% arrange(ClumpedLocus,CHR,POS) %>%
      distinct(ID,.keep_all = T)
  }
  return(out_list)
}
get_multi_nominal_loci = function(df,get_clin_nom,get_direction) {
  print("Also implementing directionality,HWE QC, and removing variants failing clumping")
  out_df = list()
  for (coh in c("NIA","NIMH","NIA_NIMH_META","UKB","AOU","UKB_AOU_META")) {
    print(paste("On cohort:",coh))
    
    # Non directionality qc
    out_df[[coh]] = df[[coh]] %>% filter(HWE_PASS,ClumpedLocus)
    if (get_clin_nom) out_df[[coh]] %<>% filter(ClinNom)

    # DIRECTIONALITY QC
    if (get_direction) {
      if (coh %in% c("NIA","NIMH","NIA_NIMH_META")) {
        rows_to_remove = numeric()
        if (nrow(out_df[[coh]]) == 0) next # NIMH failed to clump 
        for (row in 1:nrow(out_df[[coh]])) {
          curr_row = out_df[[coh]][row,]
          # varying meta signs
          if (str_detect(curr_row$NIA_NIMH_DIRECTION,"\\+") & str_detect(curr_row$NIA_NIMH_DIRECTION,"\\-")) {
            if (curr_row$NIA_P <= 0.05 & curr_row$NIMH_P <= 0.05) rows_to_remove %<>% append(row)
          }
        }
        if (length(rows_to_remove) > 0) out_df[[coh]] = out_df[[coh]][-rows_to_remove,]
      }else if (coh %in% c("UKB","AOU","UKB_AOU_META")) {
        rows_to_remove = numeric()
        for (row in 1:nrow(out_df[[coh]])) {
          curr_row = out_df[[coh]][row,]
          # varying meta signs
          if (str_detect(curr_row$UKB_AOU_DIRECTION,"\\+") & str_detect(curr_row$UKB_AOU_DIRECTION,"\\-")) {
            if (is.na(curr_row$UKB_P) | is.na(curr_row$AOU_P)) next
            if (curr_row$UKB_P <= 0.05 & curr_row$AOU_P <= 0.05) rows_to_remove %<>% append(row)
          }
        }
        if (length(rows_to_remove) > 0) out_df[[coh]] = out_df[[coh]][-rows_to_remove,]
        print(glue("Removed {length(rows_to_remove)} rows due to direction issues"))
      }
    }
  }
  return(out_df)
}
make_supplemental_table_hwe_aou_p = function(in_df) {
  midp_df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/anc_hwe_midp.txt",show_col_types=F) %>%
    mutate(ALL_GW_MAF=NA,ALL_GW_P=NA,EUR_GW_MAF=NA,EUR_GW_P = NA,
           AFR_GW_MAF=NA,AFR_GW_P = NA,AMR_GW_MAF=NA,AMR_GW_P = NA) %>% 
    separate(ID,into=c("CHR","POS","A0","A1"),sep="-",remove=F) %>%
    mutate(IDrev = glue("{CHR}-{POS}-{A1}-{A0}")) %>%
    filter(ID %in% in_df$ID | IDrev %in% in_df$ID) 
  intersect_dfs = lapply(X=c("all","eur","afr","amr"),FUN = function(x) {
    vroom(glue("single_anc_final_hits_intersect/final_hits_intersect_aou_{x}.txt"),show_col_types = F)
  })
  for (r in 1:nrow(midp_df)) {
    if (r %% 100 == 0) print(glue("On row {r} of {nrow(midp_df)}"))
    curr_r = midp_df[r,]
    for (a in 1:4) {
      match = intersect_dfs[[a]] %>% filter(ID == curr_r$ID | IDrev == curr_r$ID)
      if (nrow(match) > 1) match %<>% filter(ID == curr_r$ID)
      if (nrow(match)>0) {
        midp_df[[r,(11+(2*a))]] = match$A1FREQ
        midp_df[[r,(12+(2*a))]] = match$Pval
      }
    }
  }
  tmp = midp_df %>% rename(CHROM=CHR,Allele0=A0,Allele1=A1) %>%
    mutate(Rsid=NA,ProximalGene=NA,.after="ID") %>%
    mutate(POS = as.numeric(POS))
  out_df = getRsidsAndGenesForMissingVariants(tmp,T) %>% arrange(CHROM,POS) %>%
    select(-Allele0,-Allele1,-CHROM,-POS,-IDrev)
  return(out_df)
}
make_supplemental_table_meta_sig_not_single = function(df,df_final) {
  out_df = (df %>% filter(ALL_GW_P <= 5e-8 & EUR_GW_P > 5e-8 & AFR_GW_P > 5e-8 & 
                         AMR_GW_P > 5e-8,MIDP_EUR>=1e-15 & MIDP_AFR>=1e-15 & 
                         MIDP_AMR>=1e-15 & MIDP_EAS>=1e-15 & MIDP_SAS>=1e-15 & MIDP_MID>=1e-15)) %>%
    filter(ID %in% df_final$ID)
  return(out_df)
}
make_supplemental_table_sex = function(df) {
  out_df = data.frame(ID=df$ID,AllSexOR=df$OR,AllSexP=df$AoU_Pval,FemaleOR=NA,
                      FemaleP=NA,MaleOR=NA,MaleP=NA)
  sex_lists = lapply(X=c("female","male"),FUN=function(x) {
    vroom(glue("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_hits_{x}_stats.txt"),show_col_types = F) %>%
      mutate(Pval = 10^(-LOG10P))
  })
  for (r in 1:nrow(out_df)) {
    if (r %% 100 == 0) print(glue("On row {r} of {nrow(out_df)}"))
    curr_row = out_df[r,]
    for (s in 1:2) {
      match = sex_lists[[s]] %>% filter(ID == curr_row$ID)
      if (nrow(match) > 0) out_df[[r,(2+s*2)]] = exp(match$BETA)
      if (nrow(match)>0) out_df[[r,(3+s*2)]] = match$Pval
    }
  }
  return(out_df %>% filter(!is.na(AllSexP)))
}
make_supplemental_table_cc_approx = function(df) {
  df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ukb_meta_hits_approx_vs_cc_1_20_v_1_5_effect_aou.txt",show_col_types = F) %>%
    filter(!is.na(P_Approx_120)) %>% filter(ID %in% df$ID) %>% 
    separate(ID,into=c("CHROM","POS","Allele0","Allele1"),sep="-",remove=F) %>%
    mutate(Rsid=NA,ProximalGene=NA,POS=as.numeric(POS),.after="ID") 
  out_df = getRsidsAndGenesForMissingVariants(df,T) %>% select(-CHROM,-POS,-Allele0,-Allele1)
  print(glue("Num significant pre-downsample: {nrow(df %>% filter(P_Approx_120 <= 5e-8))} \\
             Num significant post-downsample: {nrow(df %>% filter(P_Approx_15 <= 5e-8))}"))
  return(out_df)
}
make_supplemental_table_bell_with_power = function() {
  power = write_power_df()
  df = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/bell_grs_hits_vs_aou_hwe.txt") %>%
    mutate(Rsid=NA,ProximalGene=NA,.after="ID") %>%
    rename(CHROM=CHR,Allele0=A0,Allele1=A1) %>%
    inner_join(power %>% select(-MAF),by="ID")
  out_df = getRsidsAndGenesForMissingVariants(df,T) %>%
    select(-CHROM,-POS,-Allele0,-Allele1,-ID_rev) %>%
    rename(BellOR=OR)
  
  return(out_df)
}
update_clumped_loci = function(df_in) {
  clump_data = vroom("working/clumps.txt",show_col_types = F)
  df = df_in %>% mutate(ID = glue("{CHROM}-{POS}-{Allele0}-{Allele1}"),IDrev = glue("{CHROM}-{POS}-{Allele1}-{Allele0}")) %>%
    filter(!(str_detect(Direction,"\\+") & str_detect(Direction,"\\-")))
  
  # First, add rsids to missing entries in df that somehow did not make it in
  rsids = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/FAVOR/FAVOR_Rsids.txt") %>%
    distinct(FAVOR_VCF,.keep_all = T)
  for (row in 1:nrow(df)) {
    if (!is.na(df$Rsid[[row]])) next
    if (df$ID[[row]] %in% rsids$FAVOR_VCF)
      df$Rsid[[row]] = rsids$Rsid[which(rsids$FAVOR_VCF == df$ID[[row]])]
  }
  
  # Then start the clumping
  out_df = df %>% mutate(ClumpedLocus = NA) 
  failed_clumping_loci = numeric()
  for (loc in unique(df$Locus)) {
    curr_loc = df %>% filter(Locus == loc)
    curr_loc_variants = curr_loc$ID %>% append(curr_loc$IDrev) 
    
    # get matched clump data
    clump_matches_id = clump_data %>% filter(ID == "-1") #empty
    clump_matches_sp2 = clump_data %>% filter(ID == "-1") #empty
    
    clump_matches_id = clump_data %>% filter(ID %in% curr_loc$ID | ID %in% curr_loc$IDrev)
    for (v in curr_loc_variants) clump_matches_sp2 %<>% add_row(clump_data %>% filter(str_detect(SP2,v)))
    clump_match = clump_matches_id %>% add_row(clump_matches_sp2)
    
    if (nrow(clump_match) < 1) {
      print(glue("Missing clump match locus {loc}. Lead var: {(curr_loc %>% arrange(MetaP))$ID[[1]]}"))
      if (nrow(curr_loc) == 1) out_df$ClumpedLocus[which(out_df$Locus == loc)] = TRUE
      else {
        out_df$ClumpedLocus[which(out_df$Locus == loc)] = FALSE
        failed_clumping_loci %<>% append(loc)
      }
      next
    }
    for (m in 1:nrow(clump_match)) { # loop through each unique locus that was present
      # maybe the loss of nominal hits is due to flipped alleles? It is not.
      curr_clumped_loc = curr_loc %>% filter(ID %in% clump_match$ID[[m]] |
                                               str_detect(clump_match$SP2[[m]],ID) |
                                               IDrev %in% clump_match$ID[[m]] |
                                               str_detect(clump_match$SP2[[m]],IDrev))
      target_ids = clump_match$ID[[m]] %>% append(str_split(clump_match$SP2[[m]],",")[[1]])
      new_locus = min(setdiff(1:length(unique(out_df$Locus)) + min(out_df$Locus), out_df$Locus)) # new number
      out_df$Locus[which(df$ID %in% target_ids | df$IDrev %in% target_ids)] = new_locus # new number
      out_df$ClumpedLocus[which(df$ID %in% target_ids | df$IDrev %in% target_ids)] = TRUE
    }
  }
  # Get the rows for loci that failed to clump
  failed_loci = failed_clumping_loci %>% append(out_df$Locus[which(is.na(out_df$ClumpedLocus))]) %>%
    unique()
  failed_clumps = out_df %>% filter(Locus %in% failed_loci)
  
  out_list = list(out_df,failed_clumps,out_df %>% group_by(Locus) %>% 
                    filter(AoU_Pval <= 0.05 & (UKB_Pval <= 0.05 | Niagads_Pval <= 0.05)) %>% 
                    ungroup(Locus))
  names(out_list) = c("Orig Data","Failed Clumps","Clumped Data")
  return(out_list)
}
make_manhattan_easy_label = function(dataset,df,lead_var,fig_num,hide_new_gene=F,only_common=F,
                                     hide_x = T,no_indels=F,hide_all_genes=F,maf_low = NA,
                                     no_multiallelics = F,maf_cut = NA) {
  if (dataset == "NIAGADS") {
    hwe_pass = vroom("anc_hwe_midp.txt") %>% filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,
                                                    MIDP_EAS > 1e-15,MIDP_EUR > 1e-15,
                                                    MIDP_MID > 1e-15,MIDP_SAS > 1e-15)
    gw_df = vroom("../Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_chrpos_p_1e-1.txt",show_col_types=F) %>%
      rename(CHR = chr,GENPOS=pos) %>% 
      mutate(CHR = ifelse(CHR == "X","23",CHR),CHR = as.numeric(CHR)) %>%
      mutate(ID = glue("{CHR}-{GENPOS}-{ALLELE0}-{ALLELE1}")) %>%
      mutate(IDrev = glue("{CHR}-{GENPOS}-{ALLELE1}-{ALLELE0}")) %>%
      filter(Pval > 1e-5 | ID %in% hwe_pass$ID | IDrev %in% hwe_pass$ID) %>%
      filter(CHR != 23)
  }else if (dataset == "NIA_NIMH_META") {
    hwe_pass = vroom("anc_hwe_midp.txt") %>% filter(MIDP_AFR > 1e-15,MIDP_AMR > 1e-15,
                                                    MIDP_EAS > 1e-15,MIDP_EUR > 1e-15,
                                                    MIDP_MID > 1e-15,MIDP_SAS > 1e-15)
    gw_df = vroom("final_meta_results/meta_niagads_nimh_nohwe_chrposrefalt_p_1e-1.TBL",show_col_types = F) %>%
      rename(Pval = "P-value",ID = MarkerName,GENPOS = POS) %>%
      mutate(Allele1 = toupper(Allele1),Allele2 = toupper(Allele2)) %>%
      mutate(ID = glue("{CHR}-{GENPOS}-{Allele2}-{Allele1}")) %>%
      mutate(IDrev = glue("{CHR}-{GENPOS}-{Allele1}-{Allele2}")) %>%
      filter(Pval > 1e-5 | ID %in% hwe_pass$ID | IDrev %in% hwe_pass$ID) %>%
      filter(CHR != 23)
  }else if (dataset == "AOU") {
    gw_df = vroom("../Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_p_1e-1.txt",show_col_types = F) %>%
      rename(CHR = CHROM,Freq1 = A1FREQ)
  }else if (dataset == "UKB_AOU") {
    nonmatched_gw_df = vroom("final_meta_results/meta_ukb_aou_nonmatching_chrposrefalt_p_1e-5.TBL",show_col_types = F) %>%
      rename(Pval = "P-value",ID = MarkerName,GENPOS = POS) %>%
      mutate(CHRPOS = glue("{CHR}-{GENPOS}"))
    gw_df = vroom("final_meta_results/meta_ukb_aou_chrposrefalt_p_1e-1.TBL",show_col_types = F) %>%
      rename(Pval = "P-value",ID = MarkerName,GENPOS = POS) %>%
      mutate(CHRPOS = glue("{CHR}-{GENPOS}")) %>%
      filter(CHRPOS %notin% nonmatched_gw_df) %>% # works because we remove multiallelics
      add_row(nonmatched_gw_df) %>%
      rename(ALLELE0 = Allele2,ALLELE1 = Allele1)
    
  }else if (dataset == "UKB") {
    gw_df = vroom("../Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_p_1e-1.regenie",show_col_types = F) %>%
      rename(CHR = `#CHROM`,Freq1 = A1FREQ) %>% filter(Pval <= 1,Pval > 0)
  }else if (dataset == "AOU_AFR") {
    gw_df = vroom("../Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_afr_summ_stats_AD_any_allchr_pval_chrpos_p_1e-1.txt",show_col_types = F) %>%
      rename(CHR = CHROM,Freq1 = A1FREQ)
  }else if (dataset == "AOU_AMR") {
    gw_df = vroom("../Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_amr_summ_stats_AD_any_allchr_pval_chrpos_p_1e-1.txt",show_col_types = F) %>%
      rename(CHR = CHROM,Freq1 = A1FREQ)
  }else if (dataset == "AOU_EUR") {
    gw_df = vroom("../Data_Links/AoU_GWAS/AncStrat_BT_SingleAnc_PCs/aou_eur_summ_stats_AD_any_allchr_pval_chrpos_p_1e-1.txt",show_col_types = F) %>%
      rename(CHR = CHROM,Freq1 = A1FREQ)
  }
  
  # for (v in 1:nrow(lead_var)) {
  #   match = which(gw_df$ID == lead_var$ID[[v]] | gw_df$ID == lead_var$IDrev[[v]])
  #   if (length(match)>0) gw_df[match,'IDNames'] = lead_var$Gene[[v]]
  # }

  gw_fig = gw_df %>% left_join(df,by = "ID") %>% mutate(In_Final_Hits = (ID %in% df$ID | ID %in% df$IDrev)) %>%
    mutate(CHR=CHR.x,IDNames = Gene) %>% 
    mutate(IDNames = ifelse(CHR==19 & NewOld=="Old" & POS >= 44800000 & POS <= 45200000 & !str_detect(IDNames,"APOE"),NA,IDNames)) %>%
    mutate(IDNames = str_replace_all(IDNames, "\\(dist.*?\\)", ""),
         IDNames = str_replace_all(IDNames," ",""),
         IDNames = gsub("\\([^)]*\\)", "", IDNames)) %>%
    mutate(IDNames = str_replace(IDNames,",NONE","")) %>%
    mutate(IDNames = str_replace(IDNames,"NONE,","")) %>%
    mutate(IDNames = ifelse(NewOld == "New" & !is.na(IDNames),glue("NL:{IDNames}"),IDNames)) %>%
    mutate(IDNames = ifelse(duplicated(IDNames),NA,IDNames))

  gc()
  if (hide_new_gene)
    gw_fig %<>% mutate(IDNames = ifelse(str_detect(IDNames,"NL"),"NL",IDNames))
  if (only_common)
    gw_fig %<>% filter(Freq1 >= 0.01, Freq1 <= 0.99)
  if (hide_x)
    gw_fig %<>% filter(CHR != 23)
  if (no_indels) {
    gw_fig %<>% filter(nchar(ALLELE1) == nchar(ALLELE0),nchar(ALLELE1) == 1) %>%
      remove_multiallelics()
  }
  if (hide_all_genes)
    gw_fig %<>% mutate(IDNames = NA)
  if (!is.na(maf_low))
    gw_fig %<>% filter(Freq1 >= maf_low)
  if (no_multiallelics)
    gw_fig %<>% remove_multiallelics()
  if (!is.na(maf_cut))
    gw_fig %<>% filter(Freq1 >= maf_cut & Freq1 <= (1-maf_cut))
  
  # REMOVE HWE FLAGGED RESULTS. Tested on all variants in all datasets with p <= 1e-5
  # Lead var was added as it was filtered by HWE, yet not including that was dropping some lead variants
  if (dataset != "NIAGADS" & dataset != "NIA_NIMH_META") {
    hwe = vroom("anc_hwe_midp.txt",show_col_types = F) %>% 
      filter(MIDP_AFR>1e-15,MIDP_AMR>1e-15,MIDP_EAS>1e-15,MIDP_EUR>1e-15,MIDP_MID>1e-15,MIDP_SAS>1e-15)
    gw_fig %<>% filter(Pval > 1e-5 | ID %in% hwe$ID | ID %in% hwe$IDrev | CHR == 23)
    gc()
  }
  highlight_colors = magma(4)
  highlight_colors = highlight_colors[highlight_colors != "#FCFDBFFF"]
  color_df <- data.frame(color = highlight_colors,entry=c(0:(length(highlight_colors)-1))) %>%
    filter(color != "#FCFDBFFF")
  print(ggplot(color_df, aes(x = entry, y = color)) +
    geom_tile(fill = color_df$color) +
    scale_fill_identity() +
    labs(x = ""))

  # Add lead variant labels that got cut off for unknown reasons
  if (dataset == "NIA_NIMH_META") {
    gw_fig %<>%
      mutate(IDNames = ifelse(ID == "1-165164548-G-A","RNU6-755P,LMX1A",IDNames)) %>%
      mutate(IDNames = ifelse(ID == "2-74159003-G-A","MOB1A",IDNames)) %>%
      mutate(IDNames = ifelse(ID == "5-59693193-G-A","PDE4D",IDNames)) %>%
      mutate(IDNames = ifelse(ID == "5-115784235-C-T","RNU2-49P,CDO1",IDNames)) %>%
      mutate(IDNames = ifelse(ID == "7-155194021-C-T","RN7SKP280,AC099552.1",IDNames)) %>%
      mutate(IDNames = ifelse(ID == "17-35209114-G-T","SLC35G3,AC022916.1",IDNames)) %>%
      mutate(IDNames = ifelse(IDNames == "BIN1",NA,IDNames)) %>%
      mutate(IDNames = ifelse(IDNames == "AL137789.1",NA,IDNames)) %>%
      mutate(IDNames = ifelse(IDNames == "RF00285,BCL3",NA,IDNames)) %>%
      mutate(IDNames = ifelse(IDNames == "BCL3",NA,IDNames)) %>%
      mutate(IDNames = ifelse(IDNames == "CBLC",NA,IDNames))
  }
  
  
  # All variant plot
  print("Now making plots")
  ggmanh::manhattan_plot(gw_fig,
                         pval.colname = "Pval",chr.colname = "CHR",
                         pos.colname = "GENPOS",highlight.colname = "In_Final_Hits",
                         highlight.col = highlight_colors[c(1,3)],color.by.highlight = T,
                         rescale=T,label.colname = "IDNames",
                         outfn = glue("Paper_Tables_Figures/Figure{fig_num}_all_labels.png"),
                         max.overlaps = 1000,label.font.size=2,signif = c(5e-08))
  gc()
  # # common variants
  # ggmanh::manhattan_plot(gw_fig %>% filter(Freq1 >= 0.01,Freq1 <= 0.99),
  #                        pval.colname = "Pval",chr.colname = "CHR",
  #                        pos.colname = "GENPOS",highlight.colname = "In_Final_Hits",
  #                        highlight.col = highlight_colors,color.by.highlight = T,
  #                        label.colname = "IDNames",rescale=T,
  #                        outfn = glue("Paper_Tables_Figures/Figure{fig_num}_common_0.01.png"),
  #                        max.overlaps = 1000,label.font.size=2,signif = c(5e-08))
  # 
  # # Rare variants
  # ggmanh::manhattan_plot(gw_fig %>% filter(Freq1 < 0.01 | Freq1 > 0.99),
  #                        pval.colname = "Pval",chr.colname = "CHR",
  #                        pos.colname = "GENPOS",highlight.colname = "In_Final_Hits",
  #                        highlight.col = highlight_colors,color.by.highlight = T,
  #                        label.colname = "IDNames",rescale=T,
  #                        outfn = glue("Paper_Tables_Figures/Figure{fig_num}_rare_0.01.png"),
  #                        max.overlaps = 1000,label.font.size=2,signif = c(5e-08))
  
  print(glue("Made figure. Path: Figure{fig_num}_x.pdf"))
  return(gw_fig)
}
make_supp_table_meta_hits = function() { # GW-sig hits in Meta
  df = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols_gw_sig.TBL",show_col_types = F) %>%
    rename(ID = MarkerName) %>% arrange(CHR,POS) %>%
    mutate(Rsid=NA,ProximalGene=NA,.after="ID") %>%
    rename(CHROM=CHR,Allele0=Allele1,Allele1=Allele2) %>% 
    mutate(Allele0 = toupper(Allele0),Allele1 = toupper(Allele1),IDrev = NULL)
  out_df = getRsidsAndGenesForMissingVariants(df,T)
  return(out_df)
}
make_supp_table_aou_hits = function() { # GW sig hits in AoU
  df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_chrpos_gwsig.txt",show_col_types=F) %>%
    arrange(CHROM,GENPOS) %>% filter(A1FREQ * N >= 20, (1-A1FREQ) * N >= 20) %>%
    mutate(EXTRA = NULL,Category = NULL,CHRPOS = NULL,TEST = NULL) %>% 
    mutate(Rsid = NA,ProximalGene=NA,.before = "ALLELE0") %>% 
    rename(POS=GENPOS,Allele0=ALLELE0,Allele1=ALLELE1)

  out_df = getRsidsAndGenesForMissingVariants(df,T)
  return(out_df)
}
make_supp_table_ukb_hits = function() {
  df = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_gwsig.regenie",show_col_types=F) %>%
    arrange(`#CHROM`,GENPOS) %>%
    mutate(EXTRA = NULL,Category = NULL,CHRPOS = NULL,TEST = NULL) %>%
    mutate(Rsid = NA,ProximalGene=NA,.before = "ALLELE0") %>%
    rename(POS=GENPOS,Allele0=ALLELE0,Allele1=ALLELE1,CHROM=`#CHROM`) %>%
    filter(!str_detect(ID,"\\*"),CHROM != 23)

  out_df = getRsidsAndGenesForMissingVariants(df,T)

  return(out_df)
}
make_supplemental_table_aou_anc_strat = function(df) {
  tables = lapply(c("afr","amr","eur"),function(x) {
    tmp = vroom(glue("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/AncStrat_BT/aou_AD_any_anc_{x}_gwas_hyphenid_pval_chrpos_gwsig.txt")) %>%
      filter(A1FREQ * N >= 20, (1-A1FREQ) * N >= 20) %>%
      select(-TEST,-CHISQ,-EXTRA,-FlippedFreq,-Category,-CHRPOS,-IDrev) %>%
    mutate(Rsid=NA,ProximalGene=NA,.after="ID") %>% 
      rename(POS=GENPOS,Allele0=ALLELE0,Allele1=ALLELE1) %>% 
      arrange(CHROM,POS)
    getRsidsAndGenesForMissingVariants(tmp,T)
  })
  names(tables) = c("AFR","AMR","EUR")
  return(tables)
}
write_power_df = function() {
  # Following example from Helsinski: https://www.mv.helsinki.fi/home/mjxpirin/GWAS_course/material/GWAS3.html#:~:text=Answer.,size%20also%20has%20larger%20power.
  bell_hits_pval = read_xlsx("Bellenguez_2022_SuppTables.xlsx",sheet = "Supplementary Table 5",skip=2) %>%
    select(Varianta,Chr.b,Positionc,Gened,`Minor/Major allele`,MAFe,`...13`,`Stage I + II`) %>% slice(-1) %>%
    separate(`Minor/Major allele`, into = c("MINOR", "MAJOR"), sep = "/") %>%
    separate(`Stage I + II`,into = c("OR"),sep="\\(") %>%
    rename_all(~ c("Rsid","CHR","POS","GENE","MINOR","MAJOR","MAF","Pval_Final","OR")) %>%
    select(-Rsid,-GENE) %>%
    mutate(ID=glue("{CHR}-{POS}-{MINOR}-{MAJOR}")) #POS GRCh38
  
  power = bell_hits_pval %>% select(ID,MAF,OR) %>% mutate(Power_AOU=NA,OR = as.numeric(str_trim(OR))) %>% drop_na(OR)
  for (row in 1:nrow(power)) {
    b = log(power$OR[[row]])
    n = 244838
    f = power$MAF[[row]]
    phi = 11290 / 244838
    power$Power_AOU[[row]] = pchisq(qchisq(5e-8, df = 1, lower = F), df = 1, ncp = 2*f*(1-f)*n*phi*(1-phi)*b^2, lower = F)
  }
  return(power)
}
produce_locus_zoom_data = function(df,lead) {
  for (row in 1:nrow(lead)) {
    print(glue("On row {row} of {nrow(lead)}"))
    fig_df = df %>% filter(Rsid == lead$Rsid[[row]])
    curr_chrom = fig_df$CHROM
    pos = fig_df$POS
    rsid = fig_df$Rsid
    
    # get meta window
    meta_results = vroom("aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols_p_0_01_locuszoom.TBL",show_col_types = F) %>%
      filter(CHR == curr_chrom,POS >= pos - 1e6,POS <= pos + 1e6) %>%
      mutate(Lead = MarkerName %in% lead$ID[[row]])
    
    vroom_write(meta_results,glue("locuses/{rsid}.txt"))
  }
}
isolate_phred = function(df,cutoff) {
  out_df = df %>% filter(EpiActivMax >= cutoff | EpiReprMax >= cutoff |
                           EpiTransMax >= cutoff) %>% 
    select(-CHROM,-POS,-Allele0,-Allele1)
  return(out_df)
}
produce_merged_supplemental_tables = function(df) { # take in clumped data
  lead_var = getNumberCommonRareLoci(df,maf_cut = 0.01,return_lead=T) %>%
    mutate(IDalt = glue("{CHROM}-{POS}-{Allele1}-{Allele0}"))
  
  # ST1 is just the number of cases by ancestry
  st2_aou_hits =  make_supp_table_aou_hits()
  st3_ukb_hits = make_supp_table_ukb_hits() 
  st4_approx_cc = make_supplemental_table_cc_approx(df) # approx vs exact, case control ratio
  st5_meta_hits = make_supp_table_meta_hits()
  st6_hwe_bellenguez = make_supplemental_table_bell_with_power() 
  st7_hwe_meta_hits = make_supplemental_table_hwe_aou_p(df) 
  st8_epi_phred = isolate_phred(df,cutoff=10)
  st9_gene_enrich_by_path = make_supplemental_table_cell_expression(df,clean_cell_columns = T,include_signs=T)
  st_10_meta_hwe_p = make_supplemental_table_hwe_aou_p(df) %>%
    filter(ID %in% lead_var$ID | ID %in% lead_var$IDalt)
  anc_hits = make_supplemental_table_aou_anc_strat(df)
  st_11_afr_hits = anc_hits[[1]] %>% arrange(CHROM,POS) %>% select(-CHROM,-POS,-Allele0,-Allele1)
  st_12_amr_hits = anc_hits[[2]] %>% arrange(CHROM,POS) %>% select(-CHROM,-POS,-Allele0,-Allele1)
  st_13_eur_hits = anc_hits[[3]] %>% arrange(CHROM,POS) %>% select(-CHROM,-POS,-Allele0,-Allele1)
  st_14_meta_sig_not_tested_anc = make_supplemental_table_meta_sig_not_single(make_supplemental_table_hwe_aou_p(df),df)
  
  st_list = list("SuppTable2" = st2_aou_hits,
                 "SuppTable3" = st3_ukb_hits,
                 "SuppTable4" = st4_approx_cc,
                 "SuppTable5" = st5_meta_hits,
                 "SuppTable6" = st6_hwe_bellenguez,
                 "SuppTable7" = st7_hwe_meta_hits,
                 "SuppTable8" = st8_epi_phred,
                 "SuppTable9" = st9_gene_enrich_by_path,
                 "SuppTable10" = st_10_meta_hwe_p,
                 "SuppTable11" = st_11_afr_hits,
                 "SuppTable12" = st_12_amr_hits,
                 "SuppTable13" = st_13_eur_hits,
                 "SuppTable14" = st_14_meta_sig_not_tested_anc)
  
  library(openxlsx)
  wb = createWorkbook()
  for (sheet_name in names(st_list)) {
    addWorksheet(wb, sheet_name)
    df = st_list[[sheet_name]]
    writeData(wb, sheet = sheet_name, x = df, startRow = 2, startCol = 1)
    data_cols <- 2:ncol(df)
    setColWidths(wb, sheet = sheet_name, cols = data_cols, widths = "auto")
  }
  # Save the workbook
  saveWorkbook(wb, "Paper_Tables_Figures/Supplemental_Tables.xlsx",overwrite = T)
}
getSharedHits_AoU_UKB_Meta = function() {
  aou_intersect = vroom("meta_aou_intersect_in_aou_and_ukb.txt",show_col_types = F) %>%
    mutate(ID=glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"),
           IDrev=glue("{CHROM}-{GENPOS}-{ALLELE1}-{ALLELE0}"))
  ukb_intersect = vroom("meta_ukb_intersect_in_aou_and_ukb.txt") %>%
    mutate(ID=glue("{`#CHROM`}-{GENPOS}-{ALLELE0}-{ALLELE1}"),
           IDrev=glue("{`#CHROM`}-{GENPOS}-{ALLELE1}-{ALLELE0}"))
  meta = vroom("aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols_p_0_01.TBL",
               show_col_types = F) %>% mutate(CHRPOS = glue("{CHR}-{POS}")) %>%
    filter(CHRPOS %in% aou_intersect$CHRPOS,CHRPOS %in% ukb_intersect$CHRPOS) 
  meta_id_match = meta %>%
    filter(ID %in% aou_intersect$ID | ID %in% aou_intersect$IDrev,
           ID %in% ukb_intersect$ID | ID %in% ukb_intersect$IDrev)
}
get_nominal_by_gw_sig = function(get_for_intersect,start_index=1) {
  # similar to meta intersect, get all the IDs and CHR-POS of variants significant in single studies,
  #   then intersect with other datasets, then produce a table comparing the results
  #   between studies.
  if (get_for_intersect) {
    print("This step can act up for weird file-read issues, hallucinating lack of columns")
    print("If get converted from errors, just run manually step by step")
    print("The returned output is not useful downstream")
    
    nia_gw_sig = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/NIAGADS_Personal/niagads_all_allchr_nohwe_formeta_p_gwsig.txt") %>%
      mutate(ID=glue("{chr}-{pos}-{ALLELE0}-{ALLELE1}"),
             IDrev = glue("{chr}-{pos}-{ALLELE1}-{ALLELE0}")) %>% 
      mutate(CHRPOS = glue("{chr}-{pos}")) %>%
      select(ID,IDrev,CHRPOS,Pval) 
    meta_nia_nimh_unmatched_hits = vroom("/n/home09/jwillett/true_lab_storage/00_AoU/final_meta_results/meta_nia_nimh_nonmatching_p_1e-5_chrposrefalt.TBL") %>%
      filter(`P-value` <= 5e-8) %>%
      rename(ID=MarkerName,Pval = `P-value`) %>% 
      mutate(ID = glue("{CHR}-{POS}-{toupper(Allele2)}-{toupper(Allele1)}")) %>%
      mutate(CHRPOS = glue("{CHR}-{POS}")) %>%
      mutate(IDrev = glue("{CHR}-{POS}-{toupper(Allele1)}-{toupper(Allele2)}")) %>% select(ID,IDrev,CHRPOS,Pval)
    meta_nia_nimh = vroom("final_meta_results/meta_niagads_nimh_nohwe_chrposrefalt_gw_sig.TBL") %>%
      rename(ID=MarkerName,Pval = `P-value`) %>% 
      mutate(ID = glue("{CHR}-{POS}-{toupper(Allele2)}-{toupper(Allele1)}")) %>%
      mutate(CHRPOS = glue("{CHR}-{POS}")) %>%
      mutate(IDrev = glue("{CHR}-{POS}-{toupper(Allele1)}-{toupper(Allele2)}")) %>% select(ID,IDrev,CHRPOS,Pval) %>%
      filter(ID %notin% meta_nia_nimh_unmatched_hits$ID & ID %notin% meta_nia_nimh_unmatched_hits$IDrev)
    
    ukb_gw_sig = vroom("../Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id_chrpos_gwsig.regenie",show_col_types = F) %>%
      mutate(IDrev = glue("{`#CHROM`}-{GENPOS}-{ALLELE1}-{ALLELE0}")) %>% 
      select(ID,IDrev,CHRPOS,Pval) %>% filter(Pval != 0)
    
    aou_gw_sig = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals_gwsig.txt",show_col_types = F) %>%
      mutate(IDrev = glue("{CHROM}-{GENPOS}-{ALLELE1}-{ALLELE0}")) %>% 
      mutate(CHRPOS = paste0(CHROM,"-",GENPOS)) %>%
      select(ID,IDrev,CHRPOS,Pval)
    
    meta_ukb_aou_unmatched_hits = vroom("final_meta_results/meta_ukb_aou_nonmatching_chrposrefalt_p_1e-5.TBL") %>%
      filter(`P-value` <= 5e-8) %>%
      rename(ID=MarkerName,Pval = `P-value`) %>% 
      mutate(ID = glue("{CHR}-{POS}-{toupper(Allele2)}-{toupper(Allele1)}")) %>%
      mutate(CHRPOS = glue("{CHR}-{POS}")) %>% 
      mutate(IDrev = glue("{CHR}-{POS}-{toupper(Allele1)}-{toupper(Allele2)}")) %>% select(ID,IDrev,CHRPOS,Pval)
    
    meta_ukb_aou = vroom("final_meta_results/meta_ukb_aou_chrposrefalt_gw_sig.TBL") %>%
      rename(ID=MarkerName,Pval = `P-value`) %>% 
      mutate(ID = glue("{CHR}-{POS}-{toupper(Allele2)}-{toupper(Allele1)}")) %>%
      mutate(CHRPOS = glue("{CHR}-{POS}")) %>% 
      mutate(IDrev = glue("{CHR}-{POS}-{toupper(Allele1)}-{toupper(Allele2)}")) %>% select(ID,IDrev,CHRPOS,Pval) %>%
      filter(ID %notin% meta_ukb_aou_unmatched_hits$ID & ID %notin% meta_ukb_aou_unmatched_hits$IDrev)
    
    for_intersect = nia_gw_sig %>% add_row(meta_nia_nimh) %>% add_row(meta_nia_nimh_unmatched_hits) %>%
      add_row(ukb_gw_sig) %>% add_row(aou_gw_sig) %>% add_row(meta_ukb_aou) %>%
      add_row(meta_ukb_aou_unmatched_hits)
    vroom_write(for_intersect,"working/gwas_gwsig_for_intersect.txt")
    
    just_ids_int = data.frame(ID=for_intersect$ID) %>% add_row(ID=for_intersect$IDrev)
    vroom_write(just_ids_int,"working/gwas_gwsig_for_intersect_ids.txt")
    print("Run intersect using scripts: https://github.com/juliandwillett/AoU_AlzheimersGWAS/blob/main/3.MetaAnalysis/F.IntersectByGWASSignal.sh")
    print("Intersect for local datasets (to make data organization faster)")
    return(just_ids_int)
  }else{
    print("Hardy filter not implemented, but flags present in the data frame")
    # get data for hits across studies into one table. Also filter by HWE, EXCEPT FOR X
    df = vroom("working/gwas_gwsig_for_intersect.txt",show_col_types = F) %>% select(-Pval,-CHRPOS) %>%
      mutate(NIA_P=NA,NIMH_P=NA,NIA_NIMH_META_P=NA,NIA_NIMH_DIRECTION=NA,UKB_P=NA,
             AOU_P=NA,UKB_AOU_META_P=NA,UKB_AOU_DIRECTION=NA,
             NIA_MAF=NA,NIMH_MAF=NA,NIA_NIMH_META_MAF=NA,
             UKB_MAF=NA,AOU_MAF=NA,UKB_AOU_META_MAF=NA,
             AOU_LESS20=NA,HWE_PASS=NA,AOU_AFR_P=NA,AOU_AFR_DIR=NA,
             AOU_AMR_P=NA,AOU_AMR_DIR=NA,AOU_EUR_P=NA,AOU_EUR_DIR=NA) %>% 
      filter(!str_detect(ID,"\\*")) %>%
      separate(ID,into=c("CHR","POS","NEA","EA"),remove=F) %>%
      #filter(nchar(NEA) == nchar(EA)) %>% remove_multiallelics() %>%
      distinct(ID,.keep_all = T)
    
    # Get intersect data
    nia_int = vroom("working/niagads_intersects_chrpos.txt",show_col_types = F)
    nimh_int = vroom("working/nimh_intersects_chrpos.txt",show_col_types = F)
    nia_nimh_nonmatching_meta_int = vroom("working/niagads_nimh_meta_intersects_nonmatching_chrpos.txt",show_col_types=F)
    nia_nimh_meta_int = vroom("working/niagads_nimh_meta_intersects_chrpos.txt",show_col_types = F) %>%
      mutate(ID = glue("{CHR}-{POS}-{Allele1}-{Allele2}"),
             IDrev = glue("{CHR}-{POS}-{Allele2}-{Allele1}")) %>%
      filter(ID %notin% nia_nimh_nonmatching_meta_int$MarkerName,
             IDrev %notin% nia_nimh_nonmatching_meta_int$MarkerName)
    ukb_int = vroom("working/ukb_intersects_chrpos.txt",show_col_types = F)
    aou_int = vroom("working/aou_intersects_chrpos.txt",show_col_types = F)
    ukb_aou_meta_nonmatching_int = vroom("working/ukb_aou_meta_intersects_nonmatched_chrpos.txt",show_col_types=F)
    ukb_aou_meta_int = vroom("working/ukb_aou_meta_intersects_chrpos.txt",show_col_types = F) %>%
      mutate(ID = glue("{CHR}-{POS}-{Allele1}-{Allele2}"),
             IDrev = glue("{CHR}-{POS}-{Allele2}-{Allele1}")) %>%
      filter(ID %notin% ukb_aou_meta_nonmatching_int$MarkerName,
             IDrev %notin% ukb_aou_meta_nonmatching_int$MarkerName)
    aou_afr_int = vroom("working/aou_afr_intersects_chrpos.txt")
    aou_amr_int = vroom("working/aou_amr_intersects_chrpos.txt")
    aou_eur_int = vroom("working/aou_eur_intersects_chrpos.txt")
    
    # Remove duplicate entries
    print('Removing duplicate entries where alleles swapped across studies')
    rows_to_remove = numeric()
    for (id in df$ID) {
      matches = df %>% filter(ID == id | IDrev == id)
      if (nrow(matches)>1) {
        # remove duplicates only, keeping one copy
        rows_to_remove %<>% append(which(df$ID %in% matches[2:nrow(matches),"ID"]))
      }
    }
    df = df[-rows_to_remove,]
    
    # Label HWE flagged variants
    hwe = vroom("anc_hwe_midp.txt",show_col_types = F) %>%
      filter(MIDP_EUR > 1e-15, MIDP_AFR > 1e-15, MIDP_AMR > 1e-15,
               MIDP_SAS > 1e-15, MIDP_EAS > 1e-15, MIDP_MID > 1e-15)
    df %<>% mutate(HWE_PASS = (ID %in% hwe$ID | ID %in% hwe$IDrev))

    # Fill in the table. Remove instances of failed matching
    nia_nimh_fused = nia_nimh_meta_int %>% bind_rows(nia_nimh_nonmatching_meta_int) %>%
      filter(IDrev %notin% c('19-44877713-g-t','19-44895007-c-t','19-44905579-g-t',
                             '19-44739483-g-a','19-44832419-g-a','7-152415704-c-a'))
    # skip duplicates to ensure things work. Multiallelics/indels are ultimately removed anyway.
    ukb_aou_fused = ukb_aou_meta_int %>% mutate(Matched=T) %>% bind_rows(ukb_aou_meta_nonmatching_int %>% mutate(Matched=F)) %>%
      filter(IDrev %notin% c('19-44736813-C-CT','19-44829763-TA-T','19-44844403-GA-G',
                             '19-44873636-AAC-A','19-44897776-C-CA','19-44920730-C-CA',
                             '12-21662734-AC-A','16-32421947-C-CA','17-82827811-GC-G',
                             '19-44909521-C-CT','2-111295877-ATTTTTT-A','7-443821-GGGAT-G',
                             '9-42847389-AG-A','23-49095759-C-CT','23-95602224-C-CAAA'))
    
    for (row in start_index:nrow(df)) { # check for multiple rows to deal with indels
      print(glue("On row {row} of {nrow(df)}"))
      curr_row=df[row,]
      
      # NIAGADS MATCHING
      if (curr_row$ID %in% nia_int$ID | curr_row$IDrev %in% nia_int$ID) {
        tmp = match_accounting_for_indels(nia_int,curr_row)
        if (nrow(tmp)>0) {
          df$NIA_P[[row]] = tmp$Pval
          df$NIA_MAF[[row]] = tmp$A1FREQ
        }
      }
      
      # NIMH MATCHING
      if (curr_row$ID %in% nimh_int$ID | curr_row$IDrev %in% nimh_int$ID) {
        tmp = match_accounting_for_indels(nimh_int,curr_row)
        if (nrow(tmp)>0) {
          df$NIMH_P[[row]] = tmp$Pval
          df$NIMH_MAF[[row]] = tmp$A1FREQ
        }
      }
      
      # NIA-NIMH META MATCHING
      if (curr_row$ID %in% nia_nimh_fused$MarkerName | curr_row$IDrev %in% nia_nimh_fused$MarkerName) {
        tmp = match_accounting_for_indels(nia_nimh_fused %>% mutate(ID=MarkerName),curr_row)
        if (nrow(tmp)>0) {
          df$NIA_NIMH_META_P[[row]] = tmp$`P-value`
          df$NIA_NIMH_META_MAF[[row]] = tmp$Freq1
          df$NIA_NIMH_DIRECTION[[row]] = tmp$Direction
        }
      }
      
      # UKB MATCHING
      if (curr_row$ID %in% ukb_int$ID | curr_row$IDrev %in% ukb_int$ID) {
        tmp = match_accounting_for_indels(ukb_int,curr_row)
        if (nrow(tmp)>0) {
          df$UKB_P[[row]] = tmp$Pval
          df$UKB_MAF[[row]] = tmp$A1FREQ
        }
      }
      
      # AOU MATCHING
      if (curr_row$ID %in% aou_int$ID | curr_row$IDrev %in% aou_int$ID) {
        tmp = match_accounting_for_indels(aou_int,curr_row)
        if (nrow(tmp)>0) {
          df$AOU_P[[row]] = tmp$Pval
          df$AOU_MAF[[row]] = tmp$A1FREQ
          df$AOU_LESS20[[row]] = (tmp$A1FREQ * tmp$N < 20 | (1-tmp$A1FREQ) * tmp$N < 20)
        }
      }
      
      # UKB AOU META MATCHING
      if (curr_row$ID %in% ukb_aou_fused$MarkerName | curr_row$IDrev %in% ukb_aou_fused$MarkerName) {
        tmp = match_accounting_for_indels(ukb_aou_fused %>% mutate(ID=MarkerName),curr_row) %>%
          distinct(Direction,.keep_all = T)
        if (nrow(tmp)>0) {
          df$UKB_AOU_META_P[[row]] = tmp$`P-value`
          df$UKB_AOU_META_MAF[[row]] = tmp$Freq1
          df$UKB_AOU_DIRECTION[[row]] = tmp$Direction
        }
      }
      
      # AOU AFR MATCHING
      if (curr_row$ID %in% aou_afr_int$ID | curr_row$IDrev %in% aou_afr_int$ID) {
        tmp = match_accounting_for_indels(aou_afr_int,curr_row)
        if (nrow(tmp)>0) {
          df$AOU_AFR_P[[row]] = tmp$Pval
          df$AOU_AFR_DIR[[row]] = sign(tmp$BETA)
        }
      }
      
      # AOU AMR MATCHING
      if (curr_row$ID %in% aou_amr_int$ID | curr_row$IDrev %in% aou_amr_int$ID) {
        tmp = match_accounting_for_indels(aou_amr_int,curr_row)
        if (nrow(tmp)>0) {
          df$AOU_AMR_P[[row]] = tmp$Pval
          df$AOU_AMR_DIR[[row]] = sign(tmp$BETA)
        }
      }
      
      # AOU EUR MATCHING
      if (curr_row$ID %in% aou_eur_int$ID | curr_row$IDrev %in% aou_eur_int$ID) {
        tmp = match_accounting_for_indels(aou_eur_int,curr_row)
        if (nrow(tmp)>0) {
          df$AOU_EUR_P[[row]] = tmp$Pval
          df$AOU_EUR_DIR[[row]] = sign(tmp$BETA)
        }
      }
    }
    return(df)
  }
}
match_accounting_for_indels = function(dataset,row) {
  split_id = str_split(row$ID,pattern = "-")[[1]]
  if (nchar(split_id[[3]]) == nchar(split_id[[4]])) {# SNP, so less complex 
    tmp = dataset %>% filter(ID == row$ID | ID == row$IDrev)
    
    # For matching with favor hits
    if (nrow(tmp)>1 & "variant_annovar" %in% names(dataset)) {
      tmp = tmp %>% filter(!is.na(vid))
      return(tmp)
    }
    
    # For matching in df_all_hits
    if (nrow(tmp)>1) tmp = tmp %>% filter(ID == row$ID)
    return(tmp)
  }
  else # indel, so exact match necessary
    return(dataset %>% filter(ID == row$ID))
}

make_df_all_hits_table = function(df) {
  out_df = df %>% filter(is.na(AOU_LESS20) | !AOU_LESS20) %>%
    select(-POS,-NEA,-IDrev,-AOU_LESS20) %>%
    filter(CHR != 23)
  return(out_df)
}
count_reproduced_loci = function(df,comparison) {
  total = length(unique(df$Locus))
  overlap = 0
  out_df = df %>% filter(ID == "")
  
  if (length(comparison) == 1) {
    for (loc in unique(df$Locus)) {
      loc_df = df %>% filter(Locus == loc,!!sym(comparison) <= 0.05)
      if (nrow(loc_df) > 0) {
        overlap = overlap + 1
        #out_df %<>% add_row(loc_df)
        out_df %<>% add_row(match_signs_proxy(loc_df))
        #print(match_signs_proxy(loc_df))
      }
    }
  }else{
    print("More than one comparison")
    for (loc in unique(df$Locus)) {
      loc_df = df %>% filter(Locus == loc,!!sym(comparison[[1]]) <= 0.05,!!sym(comparison[[2]]) <= 0.05)
      if (nrow(loc_df) > 0)
        overlap = overlap + 1
      #out_df %<>% add_row(match_signs_proxy(loc_df))
    }
  }
  
  #print(glue("Percentage loci overlap: {round(100 * overlap / total,2)}%"))
  print(glue("Number overlapping loci: {overlap}. Num total: {total}"))
  return(out_df)
}
make_supp_table_all_loci_variants = function(df,cohort) {
  if (cohort == 'NIA_NIMH_META')
    out_df = df %>% select(ID,Locus,CHR,Rsid,Gene,NIA_NIMH_META_MAF,EA,NIA_NIMH_DIRECTION,NIA_P,
                           NIMH_P,NIA_NIMH_META_P,UKB_AOU_META_P,NewOld) %>%
      arrange(Locus) %>% rename(MAF = NIA_NIMH_META_MAF,Direction = NIA_NIMH_DIRECTION)
  else if (cohort == "UKB_AOU_META")
    out_df = df %>% select(ID,Locus,CHR,Rsid,Gene,UKB_AOU_META_MAF,EA,UKB_AOU_DIRECTION,UKB_P,
                           AOU_P,UKB_AOU_META_P,NIA_NIMH_META_P,NewOld) %>%
      arrange(Locus) %>% rename(MAF = UKB_AOU_META_MAF,Direction = UKB_AOU_DIRECTION) %>%
      mutate(Effect = NA)
  
  # Deal with syntax
  out_df = make_all_minor_alleles(out_df %>% mutate(Effect = NA),direction_index = c(1:2))
  out_df = get_rsids_only(out_df)
  
  # Output  
  return(out_df %>% 
           mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
                  Gene = str_replace_all(Gene," ",""),
                  Gene = gsub("\\([^)]*\\)", "", Gene),
                  Gene = str_replace_all(Gene,",",", ")))
}
count_genes = function(df) {
  out_df = df %>% select(Gene) %>%
    mutate(Gene = str_replace_all(Gene, "\\(dist.*?\\)", ""),
           Gene = str_replace_all(Gene," ",""),
           Gene = gsub("\\([^)]*\\)", "", Gene),
           Gene = str_replace_all(Gene,",",", ")) %>%
    distinct(Gene,.keep_all = T)
  genes = character()
  
  for (row in 1:nrow(out_df)) {
    if (!str_detect(out_df$Gene[[row]],',')) genes %<>% append(out_df$Gene[[row]])
    else {
      spl = str_split(out_df$Gene[[row]],pattern = ',')[[1]]
      spl = str_replace_all(spl,' ','')
      genes %<>% append(spl)
    }
  }
  genes = unique(genes)
  print(glue("Num unique genes: {length(genes)}"))
}
match_signs_proxy = function(df) {
  remove_cols = numeric()
  for (row in 1:nrow(df)) {
    starting_remove_len = length(remove_cols)
    s = c(str_sub(df$NIA_NIMH_DIRECTION[[row]],1,1),
          str_sub(df$NIA_NIMH_DIRECTION[[row]],2,2),
          str_sub(df$UKB_AOU_DIRECTION[[row]],1,1),
          str_sub(df$UKB_AOU_DIRECTION[[row]],2,2))
    if (s[[1]] == "?") { # by definition here, s[[2]] has to be a + or a -
      if (s[[3]] != s[[2]] & df$UKB_P[[row]] <= 0.05  & s[[3]] != "?") remove_cols %<>% append(row)
      else if (s[[4]] != s[[2]] & df$AOU_P[[row]] <= 0.05  & s[[4]] != "?") remove_cols %<>% append(row)
    }else{ # NIA has greater power, so use that as standard
      if (s[[3]] != s[[1]] & df$UKB_P[[row]] <= 0.05 & s[[3]] != "?") remove_cols %<>% append(row)
      else if (s[[4]] != s[[1]] & df$AOU_P[[row]] <= 0.05  & s[[4]] != "?") remove_cols %<>% append(row)
    }
    if (length(remove_cols) > starting_remove_len)
      print(s)
  }
  if (length(remove_cols) == 0) return(df)
  else return(df[-remove_cols,])
}
remove_multiallelics = function(df) {
  out_df = df %>% mutate(CHRPOS = glue("{CHR}-{POS}"))
  multiallelics = vroom("working/aou_multiallelics.txt",delim="\t")
  return(out_df %>% filter(CHRPOS %notin% multiallelics$CHRPOS))
}
make_consecutive_loci = function(df) {
  out_df = df %>% arrange(Chromosome,VariantID) %>%
    mutate(Locus = as.integer(factor(Locus, levels = unique(Locus)))) %>%
    arrange(Locus)
  return(out_df)
}
get_nonmatching_ids = function() {
  # NIAGADS NIMH
  niagads_nohwe %<>% mutate(CHRPOS = glue("{chr}-{pos}"))
  nimh %<>% mutate(CHRPOS = glue("{chr}-{pos}"))
  
  niagads_nohwe %<>% filter(CHRPOS %in% nimh$CHRPOS) ; nimh %<>% filter(CHRPOS %in% niagads_nohwe$CHRPOS)
  
  niagads_nonmatching_id = niagads_nohwe %>% filter(ID %notin% nimh$ID)
  nimh_nonmatching_id = nimh %>% filter(ID %notin% niagads_nohwe$ID)
  
  nimh_convert_matching_id = nimh_nonmatching_id %>% mutate(ID = glue("{chr}-{pos}-{ALLELE1}-{ALLELE0}"))
  
  vroom_write(niagads_nonmatching_id,"working/niagads_nonmatching_id_to_nimh.txt")
  vroom_write(nimh_convert_matching_id,"working/nimh_convert_matching_id_to_niagads.txt")
  
  # UKB AOU
  ukb = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/UKB_GWAS_Data/all_variants_200k_complete_p_id.regenie")
  aou = vroom("/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_w_chr23_pvals.txt")
  
  ukb %<>% mutate(CHRPOS = glue("{`#CHROM`}-{GENPOS}"))
  aou %<>% mutate(CHRPOS = glue("{CHROM}-{GENPOS}"))
  gc()
  
  ukb %<>% filter(CHRPOS %in% aou$CHRPOS) ; aou %<>% filter(CHRPOS %in% ukb$CHRPOS)
  gc()
  
  ukb_nonmatching_id = ukb %>% filter(ID %notin% aou$ID)
  aou_nonmatching_id = aou %>% filter(ID %notin% ukb$ID)
  
  aou_convert_matching_id = aou_nonmatching_id %>% mutate(ID = glue("{CHROM}-{GENPOS}-{ALLELE1}-{ALLELE0}"))
  
  vroom_write(ukb_nonmatching_id,"working/ukb_nonmatching_id_to_aou.txt")
  vroom_write(aou_convert_matching_id,"working/aou_convert_matching_id_to_ukb.txt")
}
