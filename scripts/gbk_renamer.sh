#!/usr/bin/env bash

# path to target dir of .gbk files for clinker
dir_base=$1
# path to HydDB classification output
# NOTE: field 2 should be the HydDB group name
group_map=$2

# loop over files in dir
for gbk in "${dir_base}"/*.gbk; do
    if [[ -f $gbk ]]; then

        # get stem name and remove suffix
        gbkname="${gbk##*/}"
        gbkname="${gbkname%.*}"
        gbkname="${gbkname%-*}"

        # split the name into the genome and gene IDs
        genome_id=$(echo "$gbkname" | awk -F "___" '{print $1}')
        gene_id=$(echo "$gbkname" | awk -F "___" '{print $2}')

        # get the gene_id's corresponding HydDB group from DIAMOND output
        hydgroup=$(grep "$gene_id" "$group_map" | cut -f2)
        # species=$(grep "$genome_name" "$tax_map" | awk -F";" '{print $NF}' | sed 's/s__//' | sed 's/ /_/')

        # set new filename
        new_gbk_name="${hydgroup}___${genome_id}___${gene_id}.gbk"

        # new path name
        new_dir_name="${dir_base}/${new_gbk_name}"

        # dry run test
        echo "$new_dir_name"

        # WARN: renames files
        # mv "$gbk" "$new_dir_name"

    fi
done
