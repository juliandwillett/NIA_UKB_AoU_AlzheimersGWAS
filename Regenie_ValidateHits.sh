#########################
# Get significant hits from individual cohorts and make ID file. Run in R:
# gwas_sig = gwas %>% filter(LOG10P >= -log10(5e-8)) %>%
#   mutate(ID = glue("{CHROM}:{GENPOS}:{ALLELE0},{ALLELE1}"))
# vroom_write(gwas_sig %>% select(ID),"aou_AD_any_all_variant_anc_all_gw_hits.txt")

######################################

# get plink2
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip ;\ 
unzip plink2_linux_avx2_20221024.zip

# get regenie
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip ;\
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# get pheno/covar files. For any affected relative for AD by proxy, including grandparents (precedence for this in two papers). Use _revised file for solely 1st degree relatives.
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie_pheno.txt . ;\ 
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie_covar.txt .

# get step 1 files
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie/* .

# rename pred file as needed, ensure second string matches workspace ID.
awk '{gsub("duplicateofalzheimersgwastake5", "duplicateofalzheimersgwastake5", $2)} 1' \
    aou_step1_rg_array_anc_${anc}_pred.list > revised_pred_${anc}.list ;\


##################################

# run loop. Run by chr because faster than ancestry first. Do 1-16 and 17-22 separately due to different naming scheme

for ((i=1; i<=16; i++)); do \
  curr_chr="chr${i}" ;\

  # get curr chr files
  gsutil -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/pgen_minimal_qc/plink_${curr_chr}_multi_split* . ;\
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged \
        --geno 0.1 --mind 0.1 --hwe 1e-15 \
        --make-pgen --out plink_${curr_chr}_multi_split_merged_all \
        --set-all-var-ids @:#:\$r,\$a --new-id-max-allele-len 1015
  ./plink2 --pfile plink_${curr_chr}_multi_split_merged_all \
      --extract aou_AD_any_all_variant_anc_all_gw_hits.txt \
      --make-pgen --out plink_${curr_chr}_multi_split_merged_gw_hits
    
    # deal with loss of empty columns
    awk 'BEGIN{OFS="\t"} NR==1 {print "#FID", "IID", $2} NR>1 {print "0", $1, $2}' \
    plink_${curr_chr}_multi_split_merged_gw_hits.psam > tmp.psam ; \
    mv tmp.psam plink_${curr_chr}_multi_split_merged_gw_hits.psam
  
    # run regenie
    ./regenie_v3.2.8.gz_x86_64_Linux \
        --step 2 \
        --pgen plink_${curr_chr}_multi_split_merged_gw_hits \
        --phenoFile regenie_pheno.txt \
        --covarFile regenie_covar.txt \
        --bt \
        --firth --approx --pThresh 0.01 \
        --pred revised_pred_${anc}.list \
        --bsize 400 \
        --phenoColList AD_any \
        --out aou_step2_rg_${curr_chr}_gw_hits_anc_all \
        --minMAC 20 --mcc 
  done
  gsutil -o GSUtil:parallel_composite_upload_threshold=104857600 -m cp -r -n aou_step2_rg_* gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/rg_results_all_anc_gw_hits/ ;\
  rm plink_${curr_chr}_* ;\
done
