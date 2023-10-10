gsutil -u $GOOGLE_PROJECT -m cp -r gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/array_data_for_regenie_step1/* .

vector=("EUR" "AFR" "EAS")
for item in "${vector[@]}" ; do
    print('Curr anc:',anc)
    print(subprocess.check_output(f"./regenie_v3.2.8.gz_x86_64_Linux \
    --step 1 \
    --pgen arrays_autosomes_post_qc_pruned \
    --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar.txt \
    --bt \
    --out aou_step1_rg_array_anc_{anc} \
    --bsize 1000 \
    --keep {anc}_ids.txt \
    --lowmem \
    --lowmem-prefix tmp_rg_20", shell=True).decode('utf-8'))
    
    print(subprocess.check_output(f"gsutil -u $GOOGLE_PROJECT -m cp aou_step1_rg_array_20pcs_withsex_{anc}* {my_bucket}/data/regenie/", shell=True).decode('utf-8'))
done
