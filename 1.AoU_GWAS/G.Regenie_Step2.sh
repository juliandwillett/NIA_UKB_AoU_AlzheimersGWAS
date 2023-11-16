# Get necessary files
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip
unzip plink2_linux*
wget https://github.com/rgcgithub/regenie/releases/download/v3.2.8/regenie_v3.2.8.gz_x86_64_Linux.zip
unzip regenie_v3.2.8.gz_x86_64_Linux.zip

# Get files for current chromosome
bucket="gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67"
curr_chr="chr19"
gsutil -m cp -rn $bucket/data/pgen_minimal_qc/plink_${curr_chr}_* .

# Do QC, if it has not been run already
./plink2 --pfile plink_${curr_chr}_multi_split \
        --geno 0.1 --mind 0.1 --hwe 1e-15 \
        --make-pgen --out plink_${curr_chr}_allvar_anc_all

# revise psam file given the empty column being dropped
awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' plink_${curr_chr}_allvar_anc_all.psam > tmp
mv tmp plink_${curr_chr}_allvar_anc_all.psam

# Get files for pheno/covar/step1
gsutil -m cp -rn $bucket/data/regenie_* .
gsutil -m cp -rn $bucket/data/regenie/* .

# Rename workspace in pred file, to enable running on parallel workspaces
awk '{gsub("duplicateofalzheimersgwastake5", "duplicateofalzheimersgwastake5", $2)} 1' aou_step1_rg_array_anc_all_pred.list > revised_pred.list

# Run regenie. I recommend the "--mcc" parameter for additional QC
./regenie_v3.2.8.gz_x86_64_Linux \
    --step 2 \
    --pgen plink_${curr_chr}_allvar_anc_all \
    --phenoFile regenie_pheno_revised.txt \
    --covarFile regenie_covar_revised.txt \
    --bt --mcc --firth-se \
    --firth --approx --pThresh 0.01 \
    --pred revised_pred.list \
    --bsize 400 \
    --out aou_step2_rg_${curr_chr}_firthallvariants \
    --minMAC 20 \
    --phenoCol AD_any

# Backup results
gsutil -m cp -rn aou_step2_rg_${curr_chr}_firthallvariants* $bucket/data/rg_results_mcc/
