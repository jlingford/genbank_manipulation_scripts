#!/usr/bin/env python3
"""
BATCH reads in multiple raw KOfam scan outputs and creates a pickled dictionary of all genes and its corresponding KOfam annotation

This script does the following steps all in one go:
    1. Convert raw KOfam scan .tbl output file to a more convenient .tsv format
    2. Converts the .tsv file from step 1 into a filtered .tsv and .pkl file, containing just the best KOfam annotation per gene ID
    3. Combines all .pkl files into one big .pkl file to be used for annotating Genbank files

Input:
    - Target directory containing all KOfam scan .tbl output files

Output:
    - Dictionary .pkl file ({gene_name: description})

      Intermediate outputs:
        - KOfam annotations .tsv file (all info)
        - KOfam annotations .tsv file (filtered to best Evalue per gene)
        - KOfam dictionary .pkl file

Purpose:
    To add annotation info to Genbank file. Useful for downstream application for plotting gene neighbourhoods (clinker + gggenes). Also useful for outputting .faa files with descriptions in the header.

Prerequisites:
    - need to have ran kofam_scan on genome .faa file already

\033[1m\033[31mWARNING:\033[0m
    Runs parallel by default. Polars steps are very MEMORY hungry and have caused OOM crashes
"""
# TODO:
# - [ ] add logging
# - [x] read gzipped .tbl and .gbk files
# - [x] suppress biopython warnings
# - [x] improve output dir default
# - [x] ability to overwrite input genbank file?
# - [x] add concurrent.futures (for new script)?
# - [ ] replace pickle method with DUCKDB for better memory efficiency
# - [ ] check how the added description works with the gbk_to_faa.py script
# - [ ] better dict output from polars --> add more keys to features.qualifiers
# - [ ] option to filter out high evalue KO annotations

from io import TextIOWrapper
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
        "-i",
        "--kofam_indir",
        dest="input_dir",
        type=Path,
        metavar="DIR",
        required=True,
        help="Path to target directory containing raw KOfam scan annotation results (./*.tbl files) [Required]",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        type=Path,
        required=False,
        default="./KOfam_annotations",
        metavar="DIR",
        help="Path output directory for faa files [Default: $(pwd)/KOfam_annotations/ ]",
    )

    parser.add_argument(
        "-c",
        "--cpu",
        dest="cpu",
        type=int,
        default=None,
        metavar="N",
        required=False,
        help="No. of CPUs to use for parallelism [Default: max available]",
    )

    parser.add_argument(
        "--no_intermed",
        dest="no_intermed",
        action="store_true",
        help="If set, will not keep intermediate .tsv files from KOfam file format conversion and filtering [Default: off, will write intermediate files]",
    )

    parser.add_argument(
        "--reuse_intermed",
        dest="reuse_intermed",
        action="store_true",
        help="If intermediate files already exists, do not overwrite them. Script will read info from them instead, thus speeding up script on reruns [Default: off, writes intermediate files from scratch]",
    )

    parser.add_argument(
        "--no_parallel",
        dest="no_parallel",
        action="store_true",
        help="Will synchronously, not in parallel. Good for few inputs, very slow for many inputs",
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
            fields = line.strip().split()

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
) -> Path:
    """Converts reformatted KOfam tsv into a dictionary with relevant description info.

    KOfam description gets written in the format of:

        KO={KOFAM_ID}; Evalue={EVALUE}; Desc={DESCRIPTION};

    ---
    Args:
        input_kofam (Path): path to the reformatted KOfam table (.tsv format)

    Returns:
        output_pickle_dict (Path): path to the pickled dictionary of gene_name ID's and their corresponding KOfam annotations
    """
    ############## init output files, skip this block of code if files already exist! ############
    outtsv = Path(outdir) / "KOfam_filt_tsvs" / f"{input_kofam.stem}_filtered.tsv"
    outpkl = Path(outdir) / "KOfam_pkl_dicts" / f"{input_kofam.stem}_filtered.pkl"
    outtsv.parent.mkdir(parents=True, exist_ok=True)
    outpkl.parent.mkdir(parents=True, exist_ok=True)
    if outtsv.exists() and outpkl.exists():
        if args.reuse_intermed is True:
            return outpkl
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

    # OPTIONS: write intermediate tsv file, OR delete reformatted_kofam and don't write intermediate tsv file if args
    if args.no_intermed:
        input_kofam.unlink()
    else:
        df.write_csv(
            file=outtsv,
            include_header=True,
            separator="\t",
        )

    # write dict to pickle
    with open(outpkl, "wb") as p:
        pickle.dump(kofam_dict, p)

    return outpkl


# =============================================================
def combine_pkl_dicts(
    input_pkl_dir: list[Path], outdir: Path, args: argparse.Namespace
) -> None:
    """Combines many individual dictionaries into one big dict"""
    # init output file
    bigpkl = Path(outdir) / "KOFAM_combined_annotations.pkl"
    bigpkl.parent.mkdir(parents=True, exist_ok=True)

    # TODO: replace big pickle with DUCKDB

    # init big dict
    combined_dict = {}

    # loop over pkl dicts and combine them
    for pkl in input_pkl_dir:
        with open(pkl, "rb") as p:
            kofam_dict = pickle.load(p)
            combined_dict.update(kofam_dict)

    # WRITE FULL DICT
    with open(bigpkl, "wb") as p:
        pickle.dump(combined_dict, p)

    ######## cleanup at end ########
    if args.no_intermed:
        pkl_interm_dir = Path(args.outdir) / "KOfam_pkl_dicts"
        for file in pkl_interm_dir.glob("*.pkl"):
            file.unlink()
        pkl_interm_dir.rmdir()

    # logging
    print(f"Done. Written combined dictionary of KOfam annotations to {bigpkl}")


# =============================================================
def workflow_kofam_to_pkldict(
    input_kofam: Path,
    outdir: Path,
    args: argparse.Namespace,
) -> None:
    """Description"""

    # write reformatted_kofam tsv file
    kofam_tsv = reformat_raw_kofam(input_kofam=input_kofam, outdir=outdir, args=args)

    # extract relevant info with polars
    pkl_dict = condense_kofam_table(input_kofam=kofam_tsv, outdir=outdir, args=args)

    # logging
    print(f"Converted {input_kofam} to {pkl_dict}")


# =============================================================
def main() -> None:
    """Workflow:
    ---
    main
     └── ProcessPoolExecutor
          └── workflow_kofam_to_pkldict
               ├── reformat_raw_kofam
     │         └── condense_kofam_table
     └── combine_pkl_dicts
    """
    # parse args
    args = parse_args()

    # get all tbls in target dir
    kofam_tbls = sorted([k for k in Path(args.input_dir).glob("*.tbl*") if k.is_file()])

    ################## PARALLEL PROCESSING ######################

    if args.no_parallel:
        # run synchronously:
        for tbl in kofam_tbls:
            workflow_kofam_to_pkldict(input_kofam=tbl, outdir=args.outdir, args=args)
    else:
        # run in PARALLEL
        partial_workflow_kofam_to_pkldict = partial(
            workflow_kofam_to_pkldict,
            outdir=args.outdir,
            args=args,
        )
        with ProcessPoolExecutor(max_workers=args.cpu) as exe:
            list(exe.map(partial_workflow_kofam_to_pkldict, kofam_tbls))

    ##########################################################

    # get list of all pickle files generated above
    pkl_dicts = sorted([p for p in Path(args.outdir).rglob("*.pkl") if p.is_file()])

    # combine all smaller dicts into one big one
    combine_pkl_dicts(input_pkl_dir=pkl_dicts, outdir=args.outdir, args=args)


# =============================================================
if __name__ == "__main__":
    sys.exit(main())
