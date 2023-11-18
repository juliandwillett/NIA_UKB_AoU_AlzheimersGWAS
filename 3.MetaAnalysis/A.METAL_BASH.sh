# first ensure the ID columns are the same format in GWAS files, and create a P value column
awk -v OFS='\t' 'NR==1 {print $0 "\tPval"} NR>1 {$15 = 10^(-1 * $12); print }' aou_AD_any_anc_all_gwas.txt > aou_AD_any_anc_all_gwas_pvals.txt
awk -v OFS='\t' 'NR>1 {$3 = $1 ":" $2 ":" $4 "," $5}1' aou_AD_any_anc_all_gwas_pvals.txt > aou_AD_any_anc_all_gwas_pvals_ids.txt #necessary to make meta analysis work

awk -v OFS='\t' '{$NF = 10^(-1 * $12); print }' ukb_proxy_all_anc.txt > ukb_proxy_all_anc_with_pval.txt
awk -v OFS='\t' 'NR==1 { $13 = "Pval"}1' ukb_proxy_all_anc_with_pval.txt > ukb_proxy_all_anc_with_pvals.txt
awk -v OFS='\t' 'NR>1 {$3 = $1 ":" $2 ":" $4 "," $5}1' ukb_proxy_all_anc_with_pvals.txt > ukb_proxy_all_anc_with_pvals_ids.txt

# Run a script with the METAL parameters
salloc -p test --mem 80000 -t 0-02:00 -n 4
metal METAL_script.txt
