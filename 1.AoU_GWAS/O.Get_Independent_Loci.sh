# First, take the meta results you want, then produce single txt files for each locus
{R}
separate_loci = function(df) {
  out_list = list()
  for (loc in unique(df$Locus)) {
    out_list[[length(out_list) + 1]] = df %>% filter(Locus == loc)
  }
  for (l in 1:length(out_list)) {
    chr = out_list[[l]]$CHROM[[1]]
    vroom_write(out_list[[l]],glue("locus_hits/chr{chr}_{l}.txt"))
  }
  return(out_list)
}
{/R}

######################################
# Then bring this data into AoU (in a zip file)
