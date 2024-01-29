# Less population stratification in single ancestry data, so less concern of bias.
# Given the number of hits in the multiancestry data, check that data against single-ancestry data, 
#   removing hits that did not pass HWE in single-ancestry data.
# So all we need are the HWE p values for single-ancestry cohorts

mkdir hwe_testing ; mkdir hwe_testing/hardy_out/
ancestries=(eur afr amr sas eas mid)

# First get the hits from all analyses (GWAS1, GWAS2, Meta), produce bed file.
meta_hits = vroom("/n/holystore01/LABS/tanzi_lab/Users/jwillett/00_AoU/aou_ukb_allvar_meta_analysis_IDcolon_chrposrefalt_cols_gw_sig.TBL")
ukb_hits = vroom("/n/holystore01/LABS/tanzi_lab/Lab/UKBiobank/WGS_results/all_variants_200k_hyphen_ids_pvals_gwsig.regenie")
aou_hits = vroom("/n/holystore01/LABS/tanzi_lab/Lab/AoU/CommonPCs_NonMCC_Geno1e-1_MAC20/aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_gwsig.txt")
hits_ids = data.frame(ID=aou_hits$ID) %>% add_row(ID=meta_hits$MarkerName) %>% add_row(ID=ukb_hits$ID) %>% distinct()
hits_id = hits_id %>% separate(ID,into=c("CHR","POS","A0","A1"),sep="-")
bed_file = data.frame(CHR=hits_id_split$CHR,POSS = hits_id_split$POS,POSE = hits_id_split$POS)
vroom_write(bed_file,"all_study_hits.bed")

# Prepare bed file (if testing a single GWAS, i.e. AoU)
#awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_gwsig.txt > hwe_testing/aou_hits.bed

# Run the hardy calculations
for ((chr=1;chr<=22;chr++)); do \
  for anc in "${ancestries[@]}"; do \
    ./plink2 --pfile pgen_geno_1e-1_mac_20/chr${chr} --keep ancestries/${anc}_ids.txt \
      --missing --hardy midp --out hwe_testing/hardy_out/${chr}_${anc} --extract bed1 hwe_testing/aou_hits.bed ;\
  done \
done

# Then merge the output
for anc in "${ancestries[@]}"; do \
  head -n 1 hwe_testing/hardy_out/1_${anc}.hardy > hwe_testing/hardy_out/${anc}_hwe_stats.txt ;\
  for file in hwe_testing/hardy_out/*${anc}.hardy; do \
    tail -n +2 "$file" >> hwe_testing/hardy_out/${anc}_hwe_stats.txt ;\
  done \
done

# Rcode for producing a table for these results
aou_hits = vroom("aou_ad_any_anc_all_gwas_geno_1e-1_mac20_common_pcs_pvals_gwsig.txt",show_col_types=F) %>% 
    arrange(CHROM,GENPOS) 
anc_hwe = data.frame(ID=aou_hits$ID,MIDP_EUR=NA,MIDP_AFR=NA,MIDP_AMR=NA,MIDP_EAS=NA,MIDP_SAS=NA,MIDP_MID=NA) 

hwe_datasets = lapply(X = c("eur","afr","amr","eas","sas","mid"),FUN = function(x) {
    vroom(glue("hwe_testing/hardy_out/{x}_hwe_stats.txt"),show_col_types = F) %>% 
    distinct(ID,.keep_all = T)
})

for (row in 1:nrow(anc_hwe)) {
    curr_id = anc_hwe$ID[[row]]    
    for (anc in c(1:6)) {
        curr_anc_id = which(hwe_datasets[[anc]]$ID == curr_id)
        if (length(curr_anc_id) > 0)
            anc_hwe[[row,(1+anc)]] = hwe_datasets[[anc]]$MIDP[[curr_anc_id]]
    }
}
head(anc_hwe)
vroom_write(anc_hwe,"anc_hwe_midp.txt")

print(nrow(anc_hwe))
print(nrow(anc_hwe %>% filter(MIDP_EUR <= 1e-15 | MIDP_AFR <= 1e-15 | MIDP_AMR <= 1e-15 |
                             MIDP_EAS <= 1e-15 | MIDP_SAS <= 1e-15 | MIDP_MID <= 1e-15)))
print(nrow(anc_hwe %>% filter(MIDP_EUR <= 1e-15 | MIDP_AFR <= 1e-15 | MIDP_AMR <= 1e-15)))
