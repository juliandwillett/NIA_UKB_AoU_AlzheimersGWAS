# get ancestry predictions
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv .

# get plink2
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip ;\ 
unzip plink2_linux*

# get regenie
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip \
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# get pheno/covar files
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie_pheno.txt . ;\
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie_covar.txt .

# get step 1 files
mkdir step1_files
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie/* step1_files/

# run loop. Run by chr because faster than ancestry first. Do 1-16 and 17-22 separately due to different naming scheme
ancestries=("eur" "afr" "amr" "eas")
for ((i=1; i<=16; i++)); do   
  curr_chr="chr${i}"

  # get curr chr files
  gsutil -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/pgen_minimal_qc/plink_${curr_chr}_multi_split_merged* .
  
  for anc in "${ancestries[@]}" ; do \    
    echo "Curr chr: ${curr_chr}. Curr anc: ${anc}" ;\

    # get select variants
    ./plink2 --pfile plink_${curr_chr}_multi_split_merged \
            --export bgen-1.1 --out plink_{curr_chr}_multi_split_zeroFID_commonvariants_{anc} \
            --maf 0.01

    # run regenie step 2
    ./regenie_v3.2.8.gz_x86_64_Linux \
          --step 2 \
          --bgen plink_{curr_chr}_multi_split_zeroFID_commonvariants_{anc}.bgen \
          --sample plink_{curr_chr}_multi_split_zeroFID_commonvariants_{anc}.sample \
          --phenoFile regenie_pheno_zeroFID.txt \
          --covarFile regenie_covars_20pcs_withsex_zeroFID.txt \
          --bt \
          --keep {anc}_ids.txt \
          --firth --approx --pThresh 0.01 \
          --pred aou_step1_rg_array_20pcs_withsex_{anc}_pred.list \
          --bsize 400 \
          --out aou_step2_rg_{curr_chr}_common_20pcs_withsex_{anc} \
          --minMAC 100

    # backup files
    gsutil -m cp -r -n 
  done 
done
