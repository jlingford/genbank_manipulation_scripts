#!/usr/bin/env bash

# build three column table as input for gbk_extract_regeion_batch.py
#
# TIP: run at the base of input dir containing all genbank files
# make sure you have a list of your target genes already as a list
# have ripgrep installed

# inputs
GENELIST=$1
OUTDIR=$2
TSVNAME=$3

# use -d 1 to control max depth ripgrep searches to
# -F: fixed strings (very important for exact match)
# -w: only matches surrounded by "non-word" characters, i.e., digits, quotes, dashes (but not underscores)
# -g: glob for file match (used for extra specificity)
# $(pwd) will add the absolute path to the output
rg -d 1 -g "*.gbk" -Fw -f "$GENELIST" "$(pwd)" |
    sed 's/://' |
    sed 's#/locus_tag="##' |
    sed 's/"$//' |
    sed "s#\$#\t$OUTDIR#" |
    awk -vOFS="\t" '{$1=$1}1' >"$TSVNAME"
