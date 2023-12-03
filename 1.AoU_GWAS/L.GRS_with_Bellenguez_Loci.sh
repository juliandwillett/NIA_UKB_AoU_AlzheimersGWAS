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
vroom_write(

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

# Next run the GRS. Following: https://choishingwan.github.io/PRS-Tutorial/plink/
awk 'NR>1 {print $1}' grs_df_for_score.txt > study_snps.txt ;\
awk '{print $1,$4}' grs_df_for_score.txt > SNP.pvalue ;\
echo "0.001 0 0.001" > range_list;\
echo "0.05 0 0.05" >> range_list ;\
echo "0.1 0 0.1" >> range_list ;\
echo "0.2 0 0.2" >> range_list ;\
echo "0.3 0 0.3" >> range_list ;\
echo "0.4 0 0.4" >> range_list ;\
echo "0.5 0 0.5" >> range_list ;

./plink2 \
    --pfile tmp \
    --score grs_df_for_score.txt 1 2 3 header \
    --q-score-range range_list SNP.pvalue \
    --extract study_snps.txt \
    --out grs_chr${curr_chr}
