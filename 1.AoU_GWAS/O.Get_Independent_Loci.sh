# First, take the meta results you want, then produce single txt files for each locus
{R}
separate_loci = function(df) {
  out_list = list()
  for (loc in unique(df$Locus)) {
    curr_df = df %>% filter(Locus == loc)
    if (nrow(curr_df) > 1) out_list[[length(out_list) + 1]] = curr_df # only care about those where we need to check LD
  }
  for (l in 1:length(out_list)) {
    chr = out_list[[l]]$CHROM[[1]]
    vroom_write(out_list[[l]],glue("locus_hits/chr{chr}_{l}_numvar_{nrow(out_list[[l]])}.txt"))
  }
  return(out_list)
}
{/R}

######################################
# Then bring this data into AoU (in a zip file)
unzip locus_hits.zip
mkdir locus_hits/beds/
mkdir locus_hits/assoc_files/
mkdir locus_hits/plink1_files/
cd locus_hits
for file in *.txt ; do \
  awk '{print $5 "\t" $6 "\t" $6}' $file > beds/$file.bed ;\
  awk -F'\t' '{print $2 "\t" $14}' $file > assoc_files/$file.plink ;\
done
cd ..

# First do clumping to help reduce the number of hits
for file in locus_hits/beds/*.bed ; do \
  fname="${file##*/}" ;\
  fname="${fname%.*}" ;\
  fname="${fname%.*}" ;\
  chr="${fname%%_*}" ;\

  # First do clumping
  ./plink2 --pfile pgen_geno_1e-1_mac_20/${chr} --extract bed1 $file \
    --clump locus_hits/assoc_files/${fname}.txt.plink --clump-r2 0.01 --clump-id-field "ID" \
    --clump-p-field "MetaP" --out locus_hits/clumped/$fname
done
# merge output to make comparison easier
head -1 locus_hits/clumped/$fname.clumps > locus_hits/clumped/clumps.txt ;\
for file in locus_hits/clumped/*.clumps ; do \
  tail -n +2 "$file" >> locus_hits/clumped/clumps.txt ;\
done

################
# Then do conditional analysis, either with regenie or GCTA-COJO. Regenie would probably be better.
# For regenie conditional analysis.
for file in locus_hits/beds/*.bed ; do \
  fname="${file##*/}" ;\
  fname="${fname%.*}" ;\
  fname="${fname%.*}" ;\
  chr="${fname%%_*}" ;\

  ./plink2 --pfile pgen_geno_1e-1_mac_20/${chr} --extract bed1 $file \
    --make-pgen --out locus_hits/pgen/$fname
  awk 'NR==1 {print "#FID\tIID\tSEX"} NR>1 {print "0\t" $1 "\t" "NA"}' locus_hits/pgen/$fname.psam > t ;\
  mv t locus_hits/pgen/$fname.psam

  awk -F'\t' 'NR==1 || $14 < min {min=$14; line=$0} END {print line}' locus_hits/$fname.txt | \
    awk '{print $2}' > locus_hits/cond_data/$fname.cond

  #awk 'NR==FNR{arr[$2]; next} $3 in arr' locus_hits/$fname.txt locus_hits/pgen/$fname.pvar | \
   # awk 'NR>1 {print $3}' > locus_hits/cond_data/$fname.cond
  
  ./regenie_v3.2.8.gz_x86_64_Linux --step 2 \
    --pgen locus_hits/pgen/$fname --phenoFile regenie_pheno.txt \
    --covarFile regenie_covar_20pcs.txt --bt --firth-se --firth --approx \
    --pThresh 0.01 --pred revised_pred.list --bsize 400 --out locus_hits/cond_data_rg_out/$fname \
    --minMAC 20 --condition-list locus_hits/cond_data/$fname.cond --phenoCol AD_any ;
done
