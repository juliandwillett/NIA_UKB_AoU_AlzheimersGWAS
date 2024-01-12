# Get plink
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip
unzip plink2_linux*

# Get pruned array data for the reference
gsutil -u $GOOGLE_PROJECT -m cp -r {my_bucket}/data/array_data_pgen_files/arrays_autosomes_post_qc_pruned* .

# Get IDs for individuals from study in Hail dataset (more people have array data than srWGS)
gsutil -m cp -r -n $WORKSPACE_BUCKET/data/pgen_minimal_qc/plink_chr19_merged.psam .

# Use plink to get first 20 PCs for individuals in Hail matrix
./plink2 --pfile arrays_autosomes_post_qc_pruned_onlysrwgs --pca 20 approx \
--keep plink_chr19_merged.psam \
--out array_autosomes_postqc_pruned_onlysrwgs_pca_results_plink.eigenvec

# Backup results
gsutil -m cp -rn array_autosomes_postqc_pruned_onlysrwgs_pca_results_plink* $WORKSPACE_BUCKET/data/pc_data/

################################
# To get the rare 20 PCs
./plink2 \
  --pfile arrays_allchr \
  --max-maf 0.005 --geno 0.1 --hwe 1e-15 \
  --make-pgen --chr 1-22 \
  --out arrays_autosomes_rare_postqc \
  --indep-pairwise 100kb 1 0.1

./plink2 --pfile arrays_autosomes_post_qc --exclude arrays_allchr_post_qc.prune.out \
    --make-pgen --out arrays_autosomes_post_qc_pruned

gsutil -u $GOOGLE_PROJECT -m cp -r arrays_autosomes_post_qc* $WORKSPACE_BUCKET/data/array_data_for_regenie_step1/
