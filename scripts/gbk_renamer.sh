#!/usr/bin/env bash

# Renames .gbk files to make clinker output easier to read
#
# TODO:
# - [ ] add usage{}

# path to target dir of .gbk files for clinker
DIR_BASE=$1
# path to HydDB classification output
# NOTE: field 2 should be the HydDB group name
GROUP_MAP=$2

# loop over files in dir
for gbk in "${DIR_BASE}"/*.gbk; do
    if [[ -f $gbk ]]; then

        # assumes an input file name of 'GENOMEID___GENEID-genbank_region.gbk'

        # get stem name and remove suffix
        gbkname="${gbk##*/}"
        gbkname="${gbkname%.*}"
        gbkname="${gbkname%-*}"

        # split the name into the genome and gene IDs
        genome_id=$(echo "$gbkname" | awk -F "___" '{print $1}')
        gene_id=$(echo "$gbkname" | awk -F "___" '{print $2}')

        # get the gene_id's corresponding HydDB group from DIAMOND output
        hydgroup=$(grep "$gene_id" "$GROUP_MAP" | cut -f2)
        # species=$(grep "$genome_name" "$tax_map" | awk -F";" '{print $NF}' | sed 's/s__//' | sed 's/ /_/')

        # set new filename
        new_gbk_name="${hydgroup}___${genome_id}___${gene_id}.gbk"

        # new path name
        new_dir_name="${DIR_BASE}/${new_gbk_name}"

        # dry run test
        echo "$new_dir_name"

        # WARN: renames files
        mv "$gbk" "$new_dir_name"

    fi
done
