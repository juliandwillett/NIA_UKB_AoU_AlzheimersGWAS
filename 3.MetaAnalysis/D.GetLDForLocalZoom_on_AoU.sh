# Make SNPlist by file
for (file in list.files(path = "locuses",full.names = F)) {
    print(file)
    tmp = vroom(glue("locuses/{file}"),show_col_types = F) %>% 
        mutate(ID=glue("{CHR}-{POS}-{Allele1}-{Allele2}"),IDrev = glue("{CHR}-{POS}-{Allele2}-{Allele1}"))
    ids = data.frame(ID=tmp$ID %>% append(tmp$IDrev))
    vroom_write(ids,glue("locuses_snplist/{file}.snp"),col_names = F)
}

# produce files for common loci, noting the lead variant
chr=19 ;\
file="rs10401157.txt" ;\

for f in locuses_snplist/* ; do \
    fname=${f##*/} ;\
    fname=${fname%.*} ;\
    fname=${fname%.*} ;\
    chr=$(head -1 $f | cut -d '-' -f1) ;\
    ./plink2 --r2-phased --pfile pgen_geno_1e-1_mac_20/chr${chr} --ld-window-r2 0.0 --ld-snp-list $f \
        --out locuses_ld/$fname.ld --memory 100000 ;\
done
