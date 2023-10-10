# get ancestry predictions
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-aou-datasets-controlled/v7/wgs/short_read/snpindel/aux/ancestry/ancestry_preds.tsv .

# get regenie
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip \
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# get array data that has already been prepared (autosomes, pruned, QCd)
gsutil -u $GOOGLE_PROJECT -m cp -r gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/array_data_for_regenie_step1/* .

# get pheno/covar files
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie_pheno.txt .
gsutil -u $GOOGLE_PROJECT -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie_covar.txt .

# run regenie on individuals of each ancestry
ancestries=("eur" "afr" "amr" "eas" "sas" "mid")
for anc in "${ancestries[@]}" ; do \
    echo "Curr anc: $anc"

    # make ancestry specific ID files
    awk -v anc=$anc '$2 == anc { print $1 }' ancestry_preds.tsv > ${anc}_ids.txt

    # run regenie step 1
    ./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt \
    --out aou_step1_rg_array_anc_{anc} \
    --bsize 1000 \
    --keep {anc}_ids.txt \
    --lowmem \
    --lowmem-prefix tmp_rg", shell=True).decode('utf-8'))
    
    print(subprocess.check_output(f"gsutil -u $GOOGLE_PROJECT -m cp aou_step1_rg_array_20pcs_withsex_{anc}* {my_bucket}/data/regenie/", shell=True).decode('utf-8'))
done

gsutil -u $GOOGLE_PROJECT -m cp -r -n aou_step1_rg_array_anc_* gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/regenie/
