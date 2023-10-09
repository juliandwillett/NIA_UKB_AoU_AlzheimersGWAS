curr_chr="chr16"
gsutil -m cp -r -n gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/pgen_minimal_qc/plink_${curr_chr}_multi_split* . ;\
mkdir plink_${curr_chr}_corrected_ids ;\
for file in plink_${curr_chr}_multi_split.pgen/*.pgen ; do \
    prefix=${file%.pgen} ; prefix2=${prefix##*/} ; echo $prefix ; \        
    ./plink2 --pfile plink_${curr_chr}_multi_split.pgen/${prefix2} \
              --set-all-var-ids @:#:\$r,\$a --make-pgen --out plink_${curr_chr}_corrected_ids/${prefix2} \
              --new-id-max-allele-len 10000 ;\
    rm plink_${curr_chr}_multi_split.pgen/${prefix2}* ;\ 
done ;\

find plink_${curr_chr}_corrected_ids -type f -name "*.pgen" > pgen_files.txt ;\
sed -e 's/\.pgen//g' pgen_files.txt > pgen_files_prefix.txt ;\
FIRST_LINE=$(head -1 pgen_files_prefix.txt) ;\
tail -n +2 pgen_files_prefix.txt > pgen_files_prefix_two_on.txt ;\
./plink2 --pfile $FIRST_LINE --pmerge-list \
pgen_files_prefix_two_on.txt --make-pgen --out plink_${curr_chr}_multi_split_merged

./plink2 --pfile plink_${curr_chr}_multi_split_merged \
        --geno 0.1 --mind 0.1 --hwe 1e-15 --maf 0.01 \
        --make-pgen --out plink_${curr_chr}_multi_split_merged_common_anc_all ;\
awk 'BEGIN{OFS="\t"} NR==1 {print "#FID", "IID", $2} NR>1 {print "0", $1, $2}' \
plink_${curr_chr}_multi_split_merged_common_anc_all.psam > tmp.psam ; \
mv tmp.psam plink_${curr_chr}_multi_split_merged_common_anc_all.psam
awk '{gsub("duplicateofalzheimersgwastake5", "duplicateofduplicateofalzheimersgwastake6", $2)} 1' \
aou_step1_rg_array_anc_all_pred.list > revised_pred.list ; \
./regenie_v3.2.8.gz_x86_64_Linux \
        --step 2 \
        --pgen plink_${curr_chr}_multi_split_merged_common_anc_all \
        --phenoFile regenie_pheno.txt \
        --covarFile regenie_covar.txt \
        --bt \
        --firth --approx --pThresh 0.01 \
        --pred revised_pred.list \
        --bsize 400 \
        --out aou_step2_rg_${curr_chr}_common_anc_all \
        --minMAC 100 ;\
gsutil -o GSUtil:parallel_composite_upload_threshold=104857600 -m cp -r -n plink_${curr_chr}_multi_split_merged_common_anc_all.log gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/rg_results/ ;\
gsutil -o GSUtil:parallel_composite_upload_threshold=104857600 -m cp -r -n aou_step2_rg_${curr_chr}_common_anc_all* gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/rg_results/ ;\
gsutil -o GSUtil:parallel_composite_upload_threshold=104857600 -m cp -r -n plink_${curr_chr}_multi_split_merged.* gs://fc-secure-4029af59-df13-4d1b-b22c-2ae64cb3dc67/data/pgen_minimal_qc/
