#CHR done: 19

gsutil -m cp -rn $WORKSPACE_BUCKET/data/pgen_minimal_qc/plink_chr19_split_merged.p* .
mkdir region_files

R:
    df = vroom("aou_AD_any_anc_all_forclumping.txt") # file for individual cohort GWAS with p vals < 1e-5
    df = df %>% mutate(P = 10^(-LOG10P),SNP = glue("{CHROM}-{GENPOS}-{ALLELE0}-{ALLELE1}"))
    vroom_write(df,"aou_AD_any_anc_all_forclumping_colsnamed.txt")

    out_df = df %>% mutate(Region = NA, RegOrigin = NA)
  
    # first liberally assign regoins
    curr_reg = 1 ; reg_origin = 0
    for (row in 1:nrow(out_df)) {
      if (row %% 500 == 0) print(row)
      curr_row = out_df[row,]
      next_row = out_df[row+1,]
      
      out_df$Region[[row]] = curr_reg ; if (reg_origin == 0) reg_origin = curr_row$GENPOS
      out_df$RegOrigin[[row]] = reg_origin
      if (row == nrow(out_df)) { next }
      else if (curr_row$CHROM == next_row$CHROM) { # could be same locus
        if (next_row$GENPOS - curr_row$GENPOS > 500000) {curr_reg = curr_reg + 1 ; reg_origin = next_row$GENPOS }
      } else { curr_reg = curr_reg + 1 ; reg_origin = next_row$GENPOS }
    }

    for (reg in 1:curr_reg) {
        vroom_write(out_df %>% Region == reg,glue("region_files/region_{reg}.txt"))
    }
    

# now run plink on the file (it will calculate LD directly on the file)
awk 'BEGIN {OFS="\t"} NR==1 { print $0; next } {$3 = $1 "-" $2 "-" $4 "-" $5}1' plink_chr19_multi_split.pvar > tmp
mv tmp plink_chr19_multi_split.pvar

awk 'NR>1 {print $16}' aou_AD_any_anc_all_forclumping_colsnamed.txt > variant_ids.txt

plink2 --pfile plink_chr19_multi_split --geno 0.1 --mind 0.1 --hwe 1e-15 --make-bed --out plink_chr19_multi_split \
    --extract variant_ids.txt 

for ((i=1208;i<=1241;i++); do
plink --bfile plink_chr19_multi_split --clump region_files/region_${i}.txt \
    --clump-kb 500 --clump-r2 0.001 --clump-p1 1e-5 --clump-p2 1e-5 --out region_files/region_${i}_clumped
done
