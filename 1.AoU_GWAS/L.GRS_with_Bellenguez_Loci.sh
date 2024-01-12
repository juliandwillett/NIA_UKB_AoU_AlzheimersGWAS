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
curr_chr="1"
./plink2 --pfile plink_chr${curr_chr}_multi_split_merged \
--set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 10000 \
--missing --hardy --out grs_files/${curr_chr} \
--extract bed1 grs_hits_plink_pos_sorted.bed
