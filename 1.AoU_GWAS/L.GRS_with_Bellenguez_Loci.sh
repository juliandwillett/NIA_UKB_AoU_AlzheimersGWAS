# Get Plink2, for QC and GRS
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip; \
unzip plink2_linux_avx2_20221024.zip 

# Produce list of variants (BED) to focus on, also making file needed for Plink GRS score function
# R code
library(glue); library(tidyverse) ; library(vroom)
grs_hits = read_xlsx("Bellenguez_2022_SuppTables.xlsx",sheet = "Supplementary Table 32",skip=2)
grs_df = data.frame(ID=grs_hits$Chromosome,EA=grs_hits$`Risk allele`,Effect=grs_hits$Beta,Pval=grs_hits$`P value`) 
plink_pos_filter = data.frame(CHROM=grs_hits$Chromosome,POSSTART=grs_hits$Position,POSSTOP=grs_hits$Position)
grs_ids = data.frame(ID
  
vroom_write(plink_pos_filter,"grs_hits_plink_pos.bed")
vroom_write(grs_df,"grs_df_for_score.txt")

###########
# Start filtering/processing the AoU plink files for these hits. Do on AoU
if ((curr_chr>=17)) ; then
    ./plink2 --pfile plink_chr${curr_chr}_multi_split \
        --geno 0.1 --hwe 1e-15 --set-all-var-ids @:#:\$r,\$a \
        --make-pgen --out tmp --new-id-max-allele-len 10000 \
        --extract bed1 grs_hits_plink_pos.bed
else
    ./plink2 --pfile plink_chr${curr_chr}_multi_split_merged \
        --geno 0.1 --hwe 1e-15 --set-all-var-ids @:#:\$r,\$a \
        --make-pgen --out tmp --new-id-max-allele-len 10000 \
        --extract bed1 grs_hits_plink_pos.bed
fi

# Then get the Hardy values to show what is going on with many of the hits
for ((chr=1;chr<17;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy --out grs_files/${curr_chr} \
    --extract bed1 grs_hits_plink_pos_sorted.bed ;\
done
for ((chr=17;chr<=22;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy --out grs_files/${curr_chr} \
    --extract bed1 grs_hits_plink_pos_sorted.bed ;\
done
# merge the results for export
head -n 1 grs_files/chr4.hardy > grs_files/hardy_values.txt ;\
for file in grs_files/*.hardy; do \
    tail -n +2 "$file" >> grs_files/hardy_values.txt ;\
done

##############
# Given that some of the variants are present in fewer than 20 individuals, do the data processing in AoU workbench via R
alleles = vroom("grs_df_for_score.txt")
df = vroom("grs_files/hardy_values.txt") %>% arrange(`#CHROM`,ID) %>% 
    separate(ID,into = c("CHR", "POS"), sep = ":",remove = F) %>%
    mutate(IDrev = glue("{CHR}:{POS}:{AX},{A1}")) %>% filter(ID %in% alleles$ID | IDrev %in% alleles$ID) %>%
    distinct()
print(glue("Orig len: {nrow(alleles)}. Remaining len: {nrow(df)}"))
head(df)
nrow(df %>% filter(P <= 1e-15))
nrow(alleles)

#############
# Those are the stats pre-QC, so what about after QC, aside from HWE?
for ((chr=1;chr<17;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged --geno 0.1 \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy midp --out grs_files_postqc_except_hwe/${curr_chr} \
    --extract bed1 grs_hits_plink_pos_sorted.bed ;\
done ;\
for ((chr=17;chr<=22;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split --geno 0.1 \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy midp --out grs_files_postqc_except_hwe/${curr_chr} \
    --extract bed1 grs_hits_plink_pos_sorted.bed ;\
done
# merge the results for export
head -n 1 grs_files_postqc_except_hwe/chr4.hardy > grs_files_postqc_except_hwe/hardy_values.txt ;\
for file in grs_files_postqc_except_hwe/*.hardy; do \
    tail -n +2 "$file" >> grs_files_postqc_except_hwe/hardy_values.txt ;\
done

#####################
# Check this all for only controls
awk 'NR==1 {print} NR>1 && $4 == 0 {print $1 "\t" $2}' regenie_pheno.txt > just_controls.txt

#####################
# Check against results for study including rare and common PCs (way too many hits)
for ((chr=1;chr<17;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy midp --out commonrare_hits_hardymissing/${curr_chr} \
    --extract bed1 commonrarepcs_hits.bed ;\
done ;\
for ((chr=17;chr<=22;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy midp --out commonrare_hits_hardymissing/${curr_chr} \
    --extract bed1 commonrarepcs_hits.bed ;\
done
# merge the results for export
head -n 1 commonrare_hits_hardymissing/chr4.hardy > commonrare_hits_hardymissing/hardy_values.txt ;\
for file in commonrare_hits_hardymissing/*.hardy; do \
    tail -n +2 "$file" >> commonrare_hits_hardymissing/hardy_values.txt ;\
done
head -n 1 commonrare_hits_hardymissing/chr4.smiss > commonrare_hits_hardymissing/smiss_values.txt ;\
for file in commonrare_hits_hardymissing/*.smiss; do \
    tail -n +2 "$file" >> commonrare_hits_hardymissing/smiss_values.txt ;\
done
head -n 1 commonrare_hits_hardymissing/chr4.vmiss > commonrare_hits_hardymissing/vmiss_values.txt ;\
for file in commonrare_hits_hardymissing/*.vmiss; do \
    tail -n +2 "$file" >> commonrare_hits_hardymissing/vmiss_values.txt ;\
done

# Filter by missingness to ensure my formula works
for ((chr=1;chr<17;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged --geno 0.1 \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy midp --out commonrare_hits_hardymissing_nomissing/${curr_chr} \
    --extract bed1 commonrarepcs_hits.bed ;\
done ;\
for ((chr=17;chr<=22;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split --geno 0.1 \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy midp --out commonrare_hits_hardymissing_nomissing/${curr_chr} \
    --extract bed1 commonrarepcs_hits.bed ;\
done
head -n 1 commonrare_hits_hardymissing_nomissing/chr4.hardy > commonrare_hits_hardymissing_nomissing/hardy_values.txt ;\
for file in commonrare_hits_hardymissing_nomissing/*.hardy; do \
    tail -n +2 "$file" >> commonrare_hits_hardymissing_nomissing/hardy_values.txt ;\
done

# Code for doing MCC against the hits
for ((curr_chr=1;curr_chr<=22;curr_chr++)) ; do \
if ((curr_chr>=17)) ; then \
    ./plink2 --pfile plink_chr${curr_chr}_multi_split \
        --geno 0.1 \
        --make-pgen --out tmp \
        --extract bed1 sig_variants_for_mcc.bed ;\
else \
    ./plink2 --pfile plink_chr${curr_chr}_multi_split_merged \
        --geno 0.1 \
        --make-pgen --out tmp \
        --extract bed1 sig_variants_for_mcc.bed ;\
fi ;\
# revise psam file given the empty column being dropped
if ((curr_chr>=16)); then\
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' tmp.psam > t ;\
mv t tmp.psam ;
fi

./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen tmp \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --qt --force-qt --mcc --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred aou_step1_rg_array_common_rare_pcs_pred.list --ignore-pred \
    --bsize 400 \
    --out commonrarepcs_mcc_testing/sig_hits_qt_mcc_chr${curr_chr} \
    --minMAC 20 \
    --phenoCol AD_any
    done
    
head -n 1 commonrarepcs_mcc_testing/sig_hits_qt_mcc_chr10_AD_any.regenie > commonrarepcs_mcc_testing/mcc.txt ;\
for file in commonrarepcs_mcc_testing/*.regenie; do \
    tail -n +2 "$file" >> commonrarepcs_mcc_testing/mcc.txt ;\
done

#######################
# So adding 20 rare PCs does not help, so going back to using only the top 20 common PCs
# Now to get HWE and run regenie on only the Bell GRS hits
{r}
df_bell = vroom("SupplementalTable5_meta_vs_bellenguez_vs_niagads_comparison.txt",show_col_types = F) %>%
    separate(ID,into=c("CHR","POS","A0","A1"),sep = "-",remove = F) %>%
    select(ID,CHR,POS,A0,A1,Bell_MAF,Bell_P) %>% mutate(AoU_Multi_MAF=NA,AoU_Multi_P=NA,AoU_Multi_HWE=NA,
                                         AoU_EUR_MAF=NA,AoU_EUR_P=NA,AoU_EUR_HWE=NA,
                                         AoU_AFR_MAF=NA,AoU_AFR_P=NA,AoU_AFR_HWE=NA,
                                         AoU_AMR_MAF=NA,AoU_AMR_P=NA,AoU_AMR_HWE=NA)
head(df_bell)
# so need the MAF and P for hits (without HWE QC), so run plink hardy AND regenie on subcohorts
vroom_write(data.frame(CHR=df_bell$CHR,POSSTART=df_bell$POS,POSSTOP=df_bell$POS),"bell_hits.bed")

df_covar = vroom("regenie_covar.txt",show_col_types=F)
names(df_covar)[1:24]
vroom_write(df_covar[,c(1:24)],"regenie_covar_20commonpcs.txt")
{/r}

{bash}
mkdir bell_test ; mkdir bell_test/rg_out ; mkdir bell_test/hwe_out ; mkdir bell_test/hwe_out_anc ; mkdir bell_test/rg_out_anc
ancestries=(eur afr amr) # ancestries with adequate sample size
for ((chr=1;chr<=22;chr++)); do \
  	./plink2 --pfile plink_chr${chr}_allvar_anc_all \
        --make-pgen --out bell_test/chr${chr}_callqc \
        --extract bed1 bell_hits.bed ;\
 	awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' bell_test/chr${chr}_callqc.psam > t ;\
	mv t bell_test/chr${chr}_callqc.psam ;\
	./plink2 --pfile bell_test/chr${chr}_callqc \
    		--set-all-var-ids @-#-\$r-\$a --new-id-max-allele-len 10000 \
    		--missing --hardy midp --out bell_test/hwe_out/chr${chr} ;\
      ./regenie_v3.2.8.gz_x86_64_Linux --step 2 \
              --pgen bell_test/chr${chr}_callqc \
              --phenoFile regenie_pheno.txt \
              --covarFile regenie_covar_20commonpcs.txt \
              --bt --firth-se \
              --firth --approx --pThresh 0.01 \
              --pred aou_step1_rg_array_anc_all_pred.list \
              --bsize 400 --minMAC 20 --phenoCol AD_any \
              --out bell_test/rg_out/chr${chr}_callqc_anc_all_with_approx_all_samples ;\
	# do on anc strat data next
 	for anc in "${ancestries[@]}"; do \
                  ./plink2 --pfile bell_test/chr${chr}_callqc --keep ${anc}_ids.txt \
                          --set-all-var-ids @-#-\$r-\$a --new-id-max-allele-len 10000 \
                        --missing --hardy midp --out bell_test/hwe_out_anc/chr${chr}_${anc} ;\
                awk '{gsub("duplicateofalzheimersgwastake5d1", "duplicateofalzheimersgwastake5", $2)} 1' aou_step1_rg_array_anc_${anc}_pred.list > revised_pred.list
                ./regenie_v3.2.8.gz_x86_64_Linux --step 2 \
                         --pgen bell_test/chr${chr}_callqc \
                         --phenoFile regenie_pheno.txt \
                         --covarFile regenie_covar_20commonpcs.txt \
                         --bt --firth-se \
                         --firth --approx --pThresh 0.01 \
                         --pred revised_pred.list \
                         --bsize 400 \
                         --out bell_test/rg_out_anc/chr${chr}_callqc_anc_${anc}_with_approx_all_samples \
                         --minMAC 20 --keep ${anc}_ids.txt \
                         --phenoCol AD_any ;\
  	done
done
### organize summ stats into one file
# all anc
head -n 1 bell_test/rg_out/chr1_callqc_anc_all_with_approx_all_samples_AD_any.regenie > bell_test/rg_out/bell_hits.txt ;\
for file in bell_test/rg_out/*.regenie; do \
    tail -n +2 "$file" >> bell_test/rg_out/bell_hits.txt ;\
done
# eur
head -n 1 bell_test/rg_out_anc/chr1_callqc_anc_eur_with_approx_all_samples_AD_any.regenie > bell_test/rg_out_anc/bell_hits_eur.txt ;\
for file in bell_test/rg_out_anc/*eur_with_approx_all_samples_AD_any.regenie; do \
    tail -n +2 "$file" >> bell_test/rg_out_anc/bell_hits_eur.txt ;\
done
# afr
head -n 1 bell_test/rg_out_anc/chr1_callqc_anc_afr_with_approx_all_samples_AD_any.regenie > bell_test/rg_out_anc/bell_hits_afr.txt ;\
for file in bell_test/rg_out_anc/*afr_with_approx_all_samples_AD_any.regenie; do \
    tail -n +2 "$file" >> bell_test/rg_out_anc/bell_hits_afr.txt ;\
done
# amr
head -n 1 bell_test/rg_out_anc/chr1_callqc_anc_amr_with_approx_all_samples_AD_any.regenie > bell_test/rg_out_anc/bell_hits_amr.txt ;\
for file in bell_test/rg_out_anc/*amr_with_approx_all_samples_AD_any.regenie; do \
    tail -n +2 "$file" >> bell_test/rg_out_anc/bell_hits_amr.txt ;\
done

### organize hardy into one file
# all anc
head -n 1 bell_test/hwe_out/chr1.hardy > bell_test/hwe_out/hwe_stats.txt ;\
for file in bell_test/hwe_out/*.hardy; do \
    tail -n +2 "$file" >> bell_test/hwe_out/hwe_stats.txt ;\
done

# eur
head -n 1 bell_test/hwe_out_anc/chr1_eur.hardy > bell_test/hwe_out_anc/hwe_stats_eur.txt ;\
for file in bell_test/hwe_out_anc/*eur.hardy; do \
    tail -n +2 "$file" >> bell_test/hwe_out_anc/hwe_stats_eur.txt ;\
done
# afr
head -n 1 bell_test/hwe_out_anc/chr1_afr.hardy > bell_test/hwe_out_anc/hwe_stats_afr.txt ;\
for file in bell_test/hwe_out_anc/*afr.hardy; do \
    tail -n +2 "$file" >> bell_test/hwe_out_anc/hwe_stats_afr.txt ;\
done
# amr
head -n 1 bell_test/hwe_out_anc/chr1_amr.hardy > bell_test/hwe_out_anc/hwe_stats_amr.txt ;\
for file in bell_test/hwe_out_anc/*amr.hardy; do \
    tail -n +2 "$file" >> bell_test/hwe_out_anc/hwe_stats_amr.txt ;\
done

{/bash}

{r}
