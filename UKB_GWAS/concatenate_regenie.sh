##script to concatenate Step 2 outputs

#!/bin/bash

# Set the directory
DIR="/MW/WGS_chunks/WGS_runs/Step1_2_corrected"

# Navigate to the directory
dx cd $DIR

# Extract the header from the first file with "icd" in its name and add a '#' character
first_file=$(dx ls | grep "icd" | head -1)
dx download $first_file -o temp_first_file.txt
awk 'BEGIN {OFS="\t"} NR==1{print "#" $0}' temp_first_file.txt > concatenated_icd_files_Step1_2_corrected.txt

# Concatenate all files with "icd" in their name while skipping their headers
dx ls | grep "icd" | xargs -I {} dx download {} -o - | awk 'NR>1' >> concatenated_icd_files_Step1_2_corrected.txt

# Sort the concatenated file based on the ID column using tab as the delimiter
# The sort command first sorts by chromosome (using a custom order for X), then by position
sort -t $'\t' -k3,3 -s -o concatenated_icd_files_Step1_2_corrected.txt <(awk 'NR==1' concatenated_icd_files_Step1_2_corrected.txt) <(awk 'NR>1' concatenated_icd_files_Step1_2_corrected.txt | sort -t ':' -k1,1V -k2,2n)

# Cleanup
rm temp_first_file.txt
