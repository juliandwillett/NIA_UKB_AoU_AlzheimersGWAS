# perform QC on array data as per regenie UKBB documentation
plink2 \
  --pfile ukb_all_chrs-merge \
  --maf 0.01 --mac 100 --geno 0.1 --hwe 1e-15 \
  --mind 0.1 \
  --write-snplist --write-samples --no-id-header \
  --out qc_pass
