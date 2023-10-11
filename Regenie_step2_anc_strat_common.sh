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

# run loop
ancestries=("eur" "afr" "amr" "eas")
for ((i=1; i<=16; i++)); do   
  curr_chr="chr${i}" ;\
  for anc in "${ancestries[@]}" ; do \    
    echo "Curr chr: ${curr_chr}. Curr anc: ${anc}" ;\
  done 
done
