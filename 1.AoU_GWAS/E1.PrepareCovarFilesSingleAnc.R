# Before R script, get the PCs, through the following
ancestries=(eur amr afr) # ids in "ancestries" folder. Already filtered for being WGS, vs just array data
for anc in "${ancestries[@]}"; do \
  ./plink2 --pfile piezo2_work/array_data/arrays_autosomes_post_qc_pruned_common \
    --keep ancestries/${anc}_ids.txt --pca 20 approx \
    --out ancestries/array_pcs_${anc} ;\
done

#####
# In R, replace covar file PCs with ancestry-specific ones, then write to named file.
covar_template = vroom("regenie_covar.txt")

for (anc in c('eur','afr','amr')) {
  anc_pcs = vroom(glue("ancestries/array_pcs_{anc}.eigenvec"))
  covar_out = covar_template
  for (pc in 1:20) {
    covar_col = which(names(covar_col) == paste0("PC",pc))
    pc_table_col = which(names(anc_pcs) == paste0("PC",pc))
    covar_out[[covar_col]] = anc_pcs[[pc_table_col]]
  }
  vroom_write(covar_out,glue("regenie_covar_{anc}.txt"))
}
