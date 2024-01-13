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
    --missing --hardy --out grs_files_postqc_except_hwe/${curr_chr} \
    --extract bed1 grs_hits_plink_pos_sorted.bed ;\
done
for ((chr=17;chr<=22;chr++)) ; do \
  curr_chr="chr${chr}" ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split --geno 0.1 \
    --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
    --missing --hardy --out grs_files_postqc_except_hwe/${curr_chr} \
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
