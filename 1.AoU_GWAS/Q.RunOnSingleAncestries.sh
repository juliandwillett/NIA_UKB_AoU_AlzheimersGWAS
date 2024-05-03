################
# GET FOR ANC STRATIFIED RESULTS
ancestries=(eur amr afr) ;\

# Only autosomes for step 1, do chr 23 on step 2
for anc in "${ancestries[@]}"; do \
  ./plink2 --pfile array_data/arrays_allchr --maf 0.01 --mac 100 --geno 0.1 \
    --make-pgen --chr 1-22 --out array_data/arrays_${anc}_autosomes_post_qc --indep-pairwise 100kb 1 0.1 ;\
  ./plink2 --pfile array_data/arrays_${anc}_autosomes_post_qc --exclude array_data/arrays_${anc}_autosomes_post_qc.prune.out \
    --make-pgen --out array_data/arrays_${anc}_autosomes_post_qc_pruned ;\
  ./plink2 --pfile array_data/arrays_${anc}_autosomes_post_qc_pruned --pca 20 approx \
    --keep regenie_input/${anc}_ids.txt \
    --out array_data/array_data/arrays_${anc}_autosomes_post_qc_pruned_pcs ;\
done
