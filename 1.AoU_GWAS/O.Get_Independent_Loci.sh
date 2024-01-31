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

for file in locus_hits/beds/*.bed ; do \
  fname="${file##*/}" ;\
  fname="${fname%.*}" ;\
  fname="${fname%.*}" ;\
  chr="${fname%%_*}" ;\
  ./plink2 --pfile pgen_geno_1e-1_mac_20/${chr} --extract bed1 $file \
    --clump locus_hits/assoc_files/${fname}.txt.plink --clump-r2 0.01 --clump-id-field "ID" \
    --clump-p-field "MetaP" --out locus_hits/clumped/$fname
  
  #./plink2 --pfile pgen_geno_1e-1_mac_20/${chr} --extract bed1 $file \
  #  --r2-phased --ld-window-r2 0 --ld-window 1000 --ld-window-kb 5000 --out locus_hits/ld_out/$fname ;\
    
  # For single comparisons (to look at dprime)
  #./plink2 --pfile pgen_geno_1e-1_mac_20/chr1 --extract bed1 locus_hits/beds/chr1_1_numvar_16.txt.bed \
  #  --ld "1-9830445-T-C" "1-8774498-CAT-C" hwe-midp #--out locus_hits/ld_out/$fname ;\

done

# merge output to make comparison easier
head -1 locus_hits/clumped/$fname.clumps > locus_hits/clumped/clumps.txt ;\
for file in locus_hits/clumped/*.clumps ; do \
  tail -n +2 "$file" >> locus_hits/clumped/clumps.txt ;\
done
