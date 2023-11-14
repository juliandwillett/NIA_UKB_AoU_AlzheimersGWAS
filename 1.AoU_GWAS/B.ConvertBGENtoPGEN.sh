# Given the export parallel by shard method, there is a bit of extra work to merge bgen files to pgen, at the benefit of saving computational resources and time

# Get Plink2, for merging files
wget https://s3.amazonaws.com/plink2-assets/alpha3/plink2_linux_avx2_20221024.zip \
unzip plink2_linux* \

# Get the bgen file, and convert each shard to pgen. Get rid of bgen file to avoid maxing out memory
curr_chr="chr16"
gsutil -m cp -r -n $WORKSPACE_BUCKET/data/bgen_minimal_qc/plink_${curr_chr}_multi_split* . ;\
mkdir plink_${curr_chr}_corrected_ids ;\
for file in plink_${curr_chr}_multi_split.pgen/* ; do \
    prefix=${file} ; prefix2=${prefix##*/} ; echo $prefix ; \        
    ./plink2 --bgen plink_${curr_chr}_multi_split.bgen/${prefix2} \
             --sample plink_${curr_chr}_multi_split.sample \
              --set-all-var-ids @:#:\$r,\$a --make-pgen --out plink_${curr_chr}_corrected_ids/${prefix2} \
              --new-id-max-allele-len 10000 ;\
    rm plink_${curr_chr}_multi_split.bgen/${prefix2}* ;\ 
done ;\

# get all the pgen prefixes so they can be merged
find plink_${curr_chr}_corrected_ids -type f -name "*.pgen" > pgen_files.txt ;\
sed -e 's/\.pgen//g' pgen_files.txt > pgen_files_prefix.txt ;\
FIRST_LINE=$(head -1 pgen_files_prefix.txt) ;\
tail -n +2 pgen_files_prefix.txt > pgen_files_prefix_two_on.txt ;\

# do the merging
./plink2 --pfile $FIRST_LINE --pmerge-list \
pgen_files_prefix_two_on.txt --make-pgen --out plink_${curr_chr}_multi_split_merged

# I recommend backing up at this point, given the time investment. While backing up is additional time, it can be valuable for later.
gsutil -m cp -rn plink_${curr_chr}_multi_split_merged* $WORKSPACE_BUCKET/data/pgen_minimal_qc/
