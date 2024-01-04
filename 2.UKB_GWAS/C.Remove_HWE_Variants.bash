### HWE vals calculated separately, so need to filter the GWAS results for those with a pval <= 1e-15
# First, isolate lines of interest (those with a p-val greater than the value)
awk '$8 > 1e-15 {print $2}' /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/Interaction/sorted_hardy

awk '$10 <= 1e-15 {print}' /n/holystore01/LABS/tanzi_lab/Users/Mo/UKB/Interaction/sorted_hardy | wc -l
