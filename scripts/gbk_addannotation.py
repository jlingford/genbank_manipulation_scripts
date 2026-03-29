#!/usr/bin/env python3
"""
Adds KOfam scan annotations to Genbank file for more descriptive CDS info.

This script does the following steps all in one go:
    1. Convert raw KOfam scan .tbl output file to a more convenient .tsv format
    2. Converts the .tsv file from step 1 into a filtered .tsv and .pkl file, containing just the best KOfam annotation per gene ID
    3. Adds the KOfam annotation to an input Genbank file under the "annotations" key of its CDS features.qualifiers attribute

Input:
    - KOfam .tbl file
    - Genbank .gbk file

Output:
    - Genbank file with updated CDS annotations
    - KOfam annotations .tsv file (all info)
    - KOfam annotations .tsv file (filtered to best Evalue per gene)
    - Dictionary .pkl file ({gene_name: description})

Purpose:
    To add annotation info to Genbank file. Useful for downstream application for plotting gene neighbourhoods (clinker + gggenes). Also useful for outputting .faa files with descriptions in the header.

Prerequisites:
    - need to have ran kofam_scan on genome .faa file already

\033[1m\033[32mTip:\033[0m
    To parallelise:

        `parallel --jobs 4 --colsep '\\t' python gbk_addannotation.py -k {1} -g {2} -o {3} --reuse_intermed :::: index.tsv`

    Where index.tsv has three columns (1: path/to/kofam_tbl 2: path/to/genbank 3: path/to/outdir)
"""
# TODO:
# - [ ] add logging
# - [x] read gzipped .tbl and .gbk files
# - [x] suppress biopython warnings
# - [x] improve output dir default
# - [x] ability to overwrite input genbank file?
# - [ ] add concurrent.futures (for new script)?
# - [ ] check how the added description works with the gbk_to_faa.py script
# - [ ] better dict output from polars --> add more keys to features.qualifiers
# - [ ] option to filter out high evalue KO annotations

from io import TextIOWrapper
from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
import argparse
import gzip
import pickle
import logging
import polars as pl
import sys
import warnings


# =============================================================
# ignore biopython warnings, which just clog up the STDOUT
warnings.filterwarnings("ignore")


# =============================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-k",
        "--kofam_input",
        dest="input",
        type=Path,
        metavar="FILE",
        required=True,
        help="Path to raw KOfam scan annotation results (.tbl file) [Required]",
    )

    parser.add_argument(
        "-g",
        "--genbank_input",
        dest="genbank_input",
        type=Path,
        metavar="FILE",
        required=True,
        help="Path to Genbank file to add KOfam descriptions to [Required]",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        type=Path,
        required=False,
        default="./annotated_genbank",
        metavar="DIR",
        help="Path output directory for faa files [Default: $(pwd)/annoted_genbank/ ]",
    )

    parser.add_argument(
        "--reuse_intermed",
        dest="reuse_intermed",
        action="store_true",
        help="If intermediate files already exists, do not overwrite them. Script will read info from them instead, thus speeding up script on reruns [Default: off, writes intermediate files from scratch]",
    )

    parser.add_argument(
        "--modify_inplace",
        dest="modify_inplace",
        action="store_true",
        help="Overwrites the original/input Genbank files with the updated annotated Genbank [Default: off, writes new Genbank file]",
    )

    parser.add_argument(
        "--skip_existing",
        dest="skip_existing",
        action="store_true",
        help="If annotated genbank file has already been written, skips it without overwriting",
    )

    args = parser.parse_args()

    return args


# =============================================================
def open_gz(filepath: Path) -> TextIOWrapper:
    """Utility function: open file, even if it is gzipped"""
    if filepath.suffix == ".gz":
        return gzip.open(filepath, "rt")
    else:
        return open(filepath, "r")


# =============================================================
def reformat_raw_kofam(
    input_kofam: Path,
    outdir: Path,
    args: argparse.Namespace,
) -> Path:
    """Converts Kofam scan raw table output file into a .tsv that polars can actually read

    The default Kofam scan output is a weird tabular output with inconsistent spacing between fields. Makes it very hard to convert into .tsv with standard awk or sed, like making sure the description is a single column and not many different columns per each word. This function reads the original file, removes the unnessary asterisks, and joins the description string into a single column, and writes out a new tsv in a streaming fashion.
    Format of tsv output is:

        col1: gene_name
        col2: KOFAM ID
        col3: threshold
        col4: score
        col5: evalue
        col6: description

    ---
    Args:
        input_kofam (Path): path to the default kofam_scan output file (.tbl)
        outdir (Path): path to the output dir

    Returns:
        Path: path to the new .tsv file
    """
    ############ init output file, skip this block of code if file already exists! #########
    reformatted_kofam = (
        Path(outdir) / "KOfam_reformat_tsvs" / f"{input_kofam.stem}_reformat.tsv"
    )
    reformatted_kofam.parent.mkdir(parents=True, exist_ok=True)
    if reformatted_kofam.exists():
        if args.reuse_intermed is True:
            return reformatted_kofam
        else:
            reformatted_kofam.unlink()

    ################### convert kofam .tbl to .tsv #########################

    # concurrent streaming style reading and writing
    with (
        open_gz(filepath=input_kofam) as infile,
        open(reformatted_kofam, "w") as outfile,
    ):
        # write header line to out file first
        header = "gene_name\tKO\tthreshold\tscore\tevalue\tKO_description\n"
        outfile.write(header)

        # read infile line by line
        while line := infile.readline():
            # skip header lines
            if line.startswith("#"):
                continue

            # split line into list
            fields = line.rstrip().split()

            # remove asterix in some lines, not sure how to informative this info is...
            if fields[0] == "*":
                fields.pop(0)

            # join the 6th field and onwards into a single field
            final_field = " ".join(fields[5:])

            # create a new list of the first 5 fields plus the new description field
            new_fields = fields[0:5] + [final_field]

            # turn list into tsv line
            new_line = "\t".join(new_fields) + "\n"
            outfile.write(new_line)

    return reformatted_kofam


# =============================================================
def condense_kofam_table(
    input_kofam: Path,
    outdir: Path,
    args: argparse.Namespace,
) -> dict[str, str]:
    """Converts reformatted KOfam tsv into a dictionary with relevant description info.

    KOfam description gets written in the format of:

        KO={KOFAM_ID}; Evalue={EVALUE}; Desc={DESCRIPTION};

    ---
    Args:
        input_kofam (Path): path to the reformatted KOfam table (.tsv format)

    Returns:
        dict[str, str]: dictionary of gene_name ID's and their corresponding KOfam annotations
    """
    ############## init output files, skip this block of code if files already exist! ############
    outtsv = Path(outdir) / "KOfam_filt_tsvs" / f"{input_kofam.stem}_filtered.tsv"
    outpkl = Path(outdir) / "KOfam_pkl_dicts" / f"{input_kofam.stem}_filtered.pkl"
    outtsv.parent.mkdir(parents=True, exist_ok=True)
    outpkl.parent.mkdir(parents=True, exist_ok=True)
    if outtsv.exists() and outpkl.exists():
        if args.reuse_intermed is True:
            # load pre-made dict from pre-saved pickle file, return it
            with open(outpkl, "rb") as p:
                kofam_dict = pickle.load(p)
                return kofam_dict
        else:
            outpkl.unlink()
            outtsv.unlink()

    ################### condense kofam table ########################
    # load df
    df = pl.read_csv(
        input_kofam,
        separator="\t",
        has_header=True,
        null_values="-",
        infer_schema_length=10000,
        ignore_errors=True,
        quote_char=None,  # WARN: important!
    )
    # column_names = [
    #     "gene_name",
    #     "KO",
    #     "threshold",
    #     "score",
    #     "evalue",
    #     "KO_description",
    # ]

    # select only relevent rows for analysis
    df = df.select(
        [
            "gene_name",
            "KO",
            "evalue",
            "KO_description",
        ]
    )

    # tsv already sorted by eval (lowest first)
    df = df.unique(subset="gene_name", keep="first").sort(by="gene_name")

    # replace any weird characters in description col
    df = df.with_columns(
        pl.col("KO_description")
        .str.replace_all(">", " to ")
        .str.replace_all(";", ",")
        .str.replace_all("#", ".")
        .str.replace_all("&", " and ")
    )

    # TODO: filter out low evalue rows?

    # create col that combines KO, eval, desc info, and return only gene_name and newcol cols
    newdf = df.with_columns(
        pl.format(
            "KOFAM={}; Evalue={}; Desc={};",
            pl.col("KO"),
            pl.col("evalue"),
            pl.col("KO_description"),
        ).alias("combined_description")
    ).select(
        [
            "gene_name",
            "combined_description",
        ]
    )

    # TODO: create a bigger dict with dict[str, list[str]] for more options for annotating Genbank

    # convert this two col table into a dict with the idiomatic .iter_rows() polars way
    kofam_dict = dict(newdf.iter_rows())

    # write intermediate files
    df.write_csv(
        file=outtsv,
        include_header=True,
        separator="\t",
    )
    # write dict to pickle
    with open(outpkl, "wb") as p:
        pickle.dump(kofam_dict, p)

    return kofam_dict


# =============================================================
def update_genbank_descs(
    genbank_input: Path,
    kofam_dict: dict[str, str],
    outdir: Path,
    args: argparse.Namespace,
) -> None:
    """Add gene annotation information into a Genbank file

    ---
    Args:
        genbank_input (Path): path to the genbank file
        kofam_dict (dict[str, str]): dictionary of KOfam annotations ({gene_name: description})
        outdir (Path): path to the output directory to write results to
        args: other arguments

    Returns:
        None: writes new Genbank file with annotation information
    """
    ################ set output gbk file ########################
    if args.modify_inplace is True:
        outpath = genbank_input
    else:
        outpath = Path(outdir) / f"{genbank_input.stem}_updated.gbk"
        outpath.parent.mkdir(parents=True, exist_ok=True)
        # skip previously written file option
        if args.skip_existing and outpath.exists():
            print(f"Skipping: Genbank file already written: {outpath}")
            return None

    ############### update_genbank_descs ########################
    with open_gz(genbank_input) as infile, open(outpath, "w") as outfile:
        # read infile and modify features.qualifiers of genbank record
        for rec in SeqIO.parse(infile, "genbank"):
            for feat in rec.features:
                if feat.type != "CDS":
                    continue

                # feat.qualifiers is a dict[str, list[str]]; get gene_name from "locus_tag" key
                gene_name = feat.qualifiers.get("locus_tag", [0])[0]

                # get kofam annotation for gene
                kofam_desc = kofam_dict.get(gene_name, None)

                # features.qualifiers is a dict; can easily update it with an annotation key,val pair
                if kofam_desc is None:
                    # feat.qualifiers["annotation"] = "None"
                    feat.qualifiers.update({"annotation": "None"})
                else:
                    # feat.qualifiers["annotation"] = kofam_desc
                    feat.qualifiers.update({"annotation": kofam_desc})

            # write update genbank record to outfile
            SeqIO.write(rec, outfile, "genbank")

    # logging
    print(f"Written annotations to {genbank_input}")


# =============================================================
def main() -> None:
    # parse args
    args = parse_args()

    # write reformatted_kofam tsv file
    kofam_tsv = reformat_raw_kofam(
        input_kofam=args.input, outdir=args.outdir, args=args
    )

    # extract relevant info with polars
    kofam_dict = condense_kofam_table(
        input_kofam=kofam_tsv, outdir=args.outdir, args=args
    )

    update_genbank_descs(
        genbank_input=args.genbank_input,
        kofam_dict=kofam_dict,
        outdir=args.outdir,
        args=args,
    )


if __name__ == "__main__":
    sys.exit(main())
