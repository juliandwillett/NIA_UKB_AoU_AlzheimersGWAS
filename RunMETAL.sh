# first ensure the ID columns are the same format, and create a P value column
awk '{$NF = 10^(-1 * $12); print }' ukb_proxy_all_anc_common.txt > ukb_proxy_all_anc_common_with_pval.txt
awk 'NR==1 { $13 = "Pval"}1' ukb_proxy_all_anc_common_with_pval.txt > ukb_proxy_all_anc_common_with_pvals.txt
awk 'NR>1 {$3 = $1 ":" $2 ":" $4 "," $5}1' aou_AD_any_anc_all_common_variant_gwas_tsv_with_pvals.txt > aou_AD_any_anc_all_common_variant_gwas_tsv_with_pvals_ids.txt
awk 'NR>1 {$3 = $1 ":" $2 ":" $4 "," $5}1' ukb_proxy_all_anc_common_with_pvals.txt > ukb_proxy_all_anc_common_with_pvals_ids.txt

# Run METAL
/n/home13/dprokopenko/bin/metal
SCHEME STDERR

MARKER   ID
ALLELE   ALLELE1 ALLELE0
FREQ     A1FREQ
EFFECT   BETA
STDERR   SE
PVAL     Pval
PROCESS /n/home09/jwillett/true_lab_storage/00_AoU/GWAS_Data/aou_AD_any_anc_all_common_variant_gwas_tsv_with_pvals_ids.txt

MARKER   ID
ALLELE   ALLELE1 ALLELE0
FREQ     A1FREQ
EFFECT   BETA
STDERR   SE
PVAL     Pval
PROCESS /n/home09/jwillett/true_lab_storage/00_AoU/GWAS_Data/ukb_proxy_all_anc_common_with_pvals_ids.txt

ANALYZE HETEROGENEITY
