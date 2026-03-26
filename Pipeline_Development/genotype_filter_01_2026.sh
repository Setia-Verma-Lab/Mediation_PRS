#!/bin/bash
# plink_prefix="/path/to/plink_prefix_chr"
# sample_names="/path/to/sample_names.txt"
# snplist="/path/to/snplist.txt"
# plink_output_prefix="my_output"

# for chr in {1..22}; do 
#     plink --bfile ${plink_prefix}${chr} --keep ${sample_names} --extract ${snplist} --recode A --out "./${plink_output_prefix}${chr}_filtered"
#     #plink --bfile ${plink_prefix}${chr} --keep ${sample_names} --extract ${snplist} --recode A  --out "./${plink_output_prefix}_filtered"
# done 
echo "Plink Filtering Started..."
plink --bfile ${plink_prefix} --keep ${sample_names} --extract ${snplist} --export A  --export-allele ${effect_allele_file} --out "./${plink_output_prefix}_filtered"
echo "Done Making Raw Files ! :)"
if [ -f "./${plink_output_prefix}_filtered.raw" ] ; then
    echo "Raw file created successfully at ${plink_output_prefix}_filtered.raw"
else
    echo "Error: Raw file was not created."
    exit 1
fi