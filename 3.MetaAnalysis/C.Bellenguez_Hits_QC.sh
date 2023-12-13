################
# query AoU for Bellenguez hits to understand what is going on there (HWE)
CHR=2 ;\
POS=37304796 ;\
./plink2 --pfile plink_chr${CHR}_multi_split_merged --chr ${CHR} --from-bp ${POS} --to-bp ${POS} --geno 0.1 --hwe 1e-15 --freq --make-pgen --out tmp ;\
./plink2 --pfile tmp --freq ;\
cat *.afreq

# for ukb
./plink2 --pfile merged_chromosome_${CHR} --chr ${CHR} --from-bp ${POS} --to-bp ${POS} --geno 0.1 --hwe 1e-15 --freq --make-pgen --out tmp ;\
./plink2 --pfile tmp --freq ;\
cat *.afreq


# rerun AoU gwas to check hits
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' tmp.psam > tmp2 ;\
mv tmp2 tmp.psam ;\
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen tmp \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt --mcc --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred revised_pred.list \
    --bsize 400 \
    --out tmp_gwas \
    --phenoCol AD_any \
    --minMAC 20 ;\
cat tmp_gwas_AD_any.regenie

# query GWAS for hits
awk 'NR==1{print $0} NR>1 && $1==5 && $2 == 14724304 {print $0}' \
/n/home09/jwillett/true_lab_storage/Data_Links/AoU_GWAS/NON_MCC_GWAS/aou_AD_any_anc_all_gwas_pvals_ids_chrompos_firthse_allminoralleles.txt
