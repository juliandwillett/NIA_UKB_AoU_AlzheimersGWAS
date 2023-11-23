library(glue); library(tidyverse) ; library(vroom)
for (i in 1:22) {
    print(glue("On chr {i}"))
    gw = vroom(glue("aou_step2_rg_chr{i}_allvar_anc_all_AD_any.regenie")) %>% 
    mutate(ID = glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"))
    favor = read.table("favor_hits.txt")
    df = gw %>% filter(LOG10P >= 7.30102999566398 | ID %in% favor$V1)
    filter_df = data.frame(CHROM=df$CHROM,POSSTART=df$GENPOS,POSSTOP=df$GENPOS)
    vroom_write(filter_df,glue("sig_variants_chr{i}.txt"),col_names=F)
}

# bash code
curr_chr=22
gsutil -m cp -rn $WORKSPACE_BUCKET/data/pgen_minimal_qc/plink_chr${curr_chr}_* .
./plink2 --pfile plink_chr${curr_chr}_multi_split_merged \
        --geno 0.1 --mind 0.1 --hwe 1e-15 \
        --make-pgen --out tmp \
        --extract bed1 sig_variants_chr${curr_chr}.txt

# revise psam file given the empty column being dropped
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' tmp.psam > t ;\
mv t tmp.psam ;

awk '{gsub("duplicateofalzheimersgwastake5", "duplicateofalzheimersgwastake5", $2)} 1' aou_step1_rg_array_anc_all_pred.list > revised_pred.list

./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen tmp \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --qt --force-qt --mcc --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred revised_pred.list --ignore-pred \
    --bsize 400 \
    --out out/sig_hits_qt_mcc_chr${curr_chr} \
    --minMAC 20 \
    --phenoCol AD_any  

./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen tmp \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --qt --force-qt --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred revised_pred.list --ignore-pred \
    --bsize 400 \
    --out out/sig_hits_qt_nomcc_chr${curr_chr} \
    --minMAC 20 \
    --phenoCol AD_any ;

./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen tmp \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred revised_pred.list --ignore-pred \
    --bsize 400 \
    --out out/sig_hits_bt_chr${curr_chr} \
    --minMAC 20 \
    --phenoCol AD_any
