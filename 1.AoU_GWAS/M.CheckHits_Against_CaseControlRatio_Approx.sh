# First ensure existing hits file is present
head aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt

# Make folder to organize results
mkdir ccratio_approx_testing

# Make bed file for the hits
awk 'NR==1 {print "CHR\tPOS\tPOS" } NR>1 {print $1 "\t" $2 "\t" $2}' aou_AD_any_anc_all_gwas_pvals_ids_gwsig.txt > aou_hits.bed

# Produce code to iterate and accomplish these tasks. Run on QC'd data
for ((chr=1;chr<=22;chr++)) ; do
  curr_chr="chr${chr}"
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged \
        --make-pgen --out plink_${curr_chr}_allvar_anc_all \
        --mac 20 --geno 0.1 --hwe 1e-15 --extract bed1 aou_hits.bed
  
  ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen plink_${curr_chr}_allvar_anc_all \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred revised_pred.list \
    --bsize 400 \
    --out aou_step2_rg_${curr_chr}_firthallvariants \
    --minMAC 20 \
    --phenoCol AD_any
done


