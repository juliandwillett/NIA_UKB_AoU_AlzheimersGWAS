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
cohorts=("NIA" "NIMH" "NIA_NIMH_META" "UKB" "AOU" "UKB_AOU_META" "NIA_NIMH_UKB_AOU_META");

for coh in "${cohorts[@]}"; do \
  mkdir locus_hits/locus_hits_${coh}/beds/ ;\
  mkdir locus_hits/locus_hits_${coh}/assoc_files/ ;\
  mkdir locus_hits/locus_hits_${coh}/clumped/ ;\
  cd locus_hits/locus_hits_${coh} ;\
  for file in *.txt ; do \
    awk '{print $2 "\t" $3 "\t" $3}' $file > beds/$file.bed ;\
    # ID \
    awk -F'\t' '{print $1,$7, $8, $9, $8, $11, $12,$13,$15}' $file > assoc_files/$file.plink ;\
    # IDrev \
    awk -F'\t' 'NR > 1 {print $6,$7, $8, $9, $8, $11, $12,$13,$15}' $file >> assoc_files/$file.plink ;\
  done ;\
  cd ../.. ;\
done
  
# First do clumping to help reduce the number of hits
cd locus_hits
for coh in "${cohorts[@]}"; do \
  for file in locus_hits_${coh}/beds/*.bed ; do \
    fname="${file##*/}" ;\
    fname="${fname%.*}" ;\
    fname="${fname%.*}" ;\
    chr="${fname%%_*}" ;\
  
    # First do clumping
    if [[ $coh == "NIA" ]]; then\
      p_field="NIA_P" ;\
    elif [[ $coh == "NIMH" ]]; then\
      p_field="NIMH_P";\
    elif [[ $coh == "NIA_NIMH_META" ]]; then\
      p_field="NIA_NIMH_META_P" ;\
    elif [[ $coh == "UKB" ]]; then\
      p_field="UKB_P" ;\
    elif [[ $coh == "AOU" ]]; then\
      p_field="AOU_P" ;\
    elif [[ $coh == "UKB_AOU_META" ]]; then\
      p_field="UKB_AOU_META_P" ;\
    elif [[ $coh == "NIA_NIMH_UKB_AOU_META" ]]; then\
      p_field="NIA_NIMH_UKB_AOU_META_P" ;\
    fi ;\
      
    ../plink2 --pfile ../pgen_qc/${chr}_geno_mac --extract bed1 $file \
      --clump locus_hits_${coh}/assoc_files/${fname}.txt.plink --clump-r2 0.01 --clump-id-field "ID" \
      --clump-p-field $p_field --out locus_hits_${coh}/clumped/$fname ;\
  done \
done
  
# merge output to make comparison easier
for coh in "${cohorts[@]}"; do \
  for file in locus_hits_${coh}/clumped/*.clumps ; do \
    fname="${file##*/}" ;\
    fname="${fname%.*}" ;\
    fname="${fname%.*}" ;\
  done ;\
  
  if [[ $coh == "NIMH" ]]; then \
    continue ;\
  fi ;\
  head -1 locus_hits_${coh}/clumped/$fname.clumps > locus_hits_${coh}/clumped/clumps.txt ;\
  for file in locus_hits_${coh}/clumped/*.clumps ; do \
    tail -n +2 "$file" >> locus_hits_${coh}/clumped/clumps.txt ;\
  done;\
done

mkdir clumps
for coh in "${cohorts[@]}"; do \
  mv locus_hits_${coh}/clumped/clumps.txt clumps/clumps_${coh}.txt ;\
done

################


DEPRECATED CODE BELOW

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
