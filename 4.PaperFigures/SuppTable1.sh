# Code for obtaining numbers for Supplemental Table 1

ancestries=(all eur afr amr eas sas mid) ;\
echo "Counts" > supp_table_1_counts.txt ;\
for anc in "${ancestries[@]}"; do # first file is pheno file, second is IDs \ 
  echo "Iterating through anc ${anc}" ;\
  echo "CURR ANC: ${anc}. Numbers reported as TOTAL, CASES, CONTROLS" >> supp_table_1_counts.txt ;\
  awk 'NR==FNR{arr[$2]; next} $2 in arr' ancestries/${anc}_ids.txt regenie_pheno.txt | wc -l >> supp_table_1_counts.txt ;\
  awk 'NR==FNR{arr[$2]; next} $2 in arr' ancestries/${anc}_ids.txt regenie_pheno.txt | awk '$4 == 1 {print}' | wc -l >> supp_table_1_counts.txt ;\
  awk 'NR==FNR{arr[$2]; next} $2 in arr' ancestries/${anc}_ids.txt regenie_pheno.txt | awk '$4 == 0 {print}' | wc -l >> supp_table_1_counts.txt ;\
done
