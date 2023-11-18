# First, you download the array files as bgen, and then convert to pgen (to make things faster). I am starting here from pgen

gsutil -u $GOOGLE_PROJECT -m cp -r $WORKSPACE_BUCKET/data/array_data_pgen_files/arrays_allchr* .

# Do QC
./plink2 \
  --pfile arrays_allchr \
  --maf 0.01 --mac 100 --geno 0.01 --hwe 1e-15 \
  --mind 0.05 \ 
  --make-pgen --chr 1-22 \
  --out arrays_autosomes_post_qc \
  --indep-pairwise 100kb 1 0.1
  
  ./plink2 --pfile arrays_autosomes_post_qc --exclude arrays_allchr_post_qc.prune.out \
    --make-pgen --out arrays_autosomes_post_qc_pruned

# Backup results
gsutil -u $GOOGLE_PROJECT -m cp -r arrays_autosomes_post_qc* $WORKSPACE_BUCKET/data/array_data_for_regenie_step1/
