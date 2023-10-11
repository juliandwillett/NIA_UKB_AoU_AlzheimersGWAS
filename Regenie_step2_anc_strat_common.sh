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
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie/* .

# rename pred file as needed, ensure second string matches workspace ID.
ancestries=("eur" "afr" "amr") ;\
for anc in "${ancestries[@]}" ; do \ 
  awk '{gsub("duplicateofalzheimersgwastake5d1", "duplicateofalzheimersgwastake5", $2)} 1' \
      step1_files/aou_step1_rg_array_anc_${anc}_pred.list > step1_files/revised_pred_${anc}.list ;\
done


##################################

# run loop. Run by chr because faster than ancestry first. Do 1-16 and 17-22 separately due to different naming scheme

for ((i=1; i<=16; i++)); do \
  curr_chr="chr${i}" ;\

  # get curr chr files
  gsutil -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/pgen_minimal_qc/plink_${curr_chr}_multi_split_merged* . ;\
  
  for anc in "${ancestries[@]}" ; do \    
    echo "Curr chr: ${curr_chr}. Curr anc: ${anc}" ;\
    awk -v anc=$anc '$2 == anc { print "0" "\t" $1 }' ancestry_preds.tsv > ${anc}_ids.txt ;\

    # run
    ./plink2 --pfile plink_${curr_chr}_multi_split_merged \
        --geno 0.1 --mind 0.1 --hwe 1e-15 --maf 0.01 \
        --make-pgen --out plink_${curr_chr}_multi_split_merged_common_anc_${anc} ;\

    # deal with loss of empty columns
    awk 'BEGIN{OFS="\t"} NR==1 {print "#FID", "IID", $2} NR>1 {print "0", $1, $2}' \
    plink_${curr_chr}_multi_split_merged_common_anc_${anc}.psam > tmp.psam ; \
    mv tmp.psam plink_${curr_chr}_multi_split_merged_common_anc_${anc}.psam ;\

    # run regenie
    ./regenie_v3.2.8.gz_x86_64_Linux \
        --step 2 \
        --pgen plink_${curr_chr}_multi_split_merged_common_anc_${anc} \
        --phenoFile regenie_pheno.txt \
        --covarFile regenie_covar.txt \
        --bt \
        --keep ${anc}_ids.txt \
        --firth --approx --pThresh 0.01 \
        --pred step1_files/revised_pred_${anc}.list \
        --bsize 400 \
        --out aou_step2_rg_${curr_chr}_common_anc_${anc} \
        --minMAC 100
  done 
done

# backup results
gsutil -o GSUtil:parallel_composite_upload_threshold=104857600 -m cp -r -n aou_step2_rg_${curr_chr}_common_anc_all* gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/rg_results/ ;\
