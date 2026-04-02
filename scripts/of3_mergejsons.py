#!/usr/bin/env python3
"""
Merge many different .json query files into one big query json for batch OpenFold3

Input:
    - protein fasta file (.faa | .fasta)
        - assumes that one fasta file = one OF3 query input. i.e., a fasta file of multiple sequences is a multimer of different chains
        - will not work with the ColabFold multimer format (i.e., multimer chains separated by ":")
Output:
    - OF3 query .json file
        - one query only; query named after the input fasta filename
        - .json file named after the input .faa file

\033[1m\033[31mWARNING:\033[0m
    ...

"""
# TODO:
# - [ ] add gzip .faa file reading?
# - [ ] improve chain labels to be infinite?
# - [ ] parse for identical sequences and make chain_ids a list[str]
# - [ ] parse description for msa.a3m filepath info and add to msa field of dict

from pathlib import Path
from typing import TextIO, NamedTuple
import argparse
import gzip
import sys
import json
from datetime import datetime


# =============================================================================
# CLI args
# =============================================================================
class Args(NamedTuple):
    indir: Path
    outpath: Path


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--indir",
        type=Path,
        metavar="FILE",
        required=True,
        help="Path to input [Required]",
    )

    parser.add_argument(
        "-o",
        "--outpath",
        type=Path,
        required=False,
        default="./OF3_merged_queries.json",
        metavar="DIR",
        help="Path/to/output_file [Optional][Default: OF3_merged_queries.json]",
    )

    args = Args(**vars(parser.parse_args()))

    return args


# =============================================================================
# Core func.
# =============================================================================
def merge_many_jsons(
    json_list: list[Path],
    outpath: Path,
    args: Args,
) -> None:
    """Converts protein fasta info to json format for OpenFold3

    ---
    Args:
        input_faa (Path): path to .faa file
        outdir (Path): path to output dir
        args (Args): other args

    Returns:
        None: Writes out .json file
    """
    ########### set outfile #############
    datestring = datetime.now().strftime("%Y-%m-%d_%H-%M-%S")
    outpath = Path(outpath.parent) / f"{outpath.stem}-{datestring}.json"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    ########### parse json files ###########

    merged_queries = {}

    for j in json_list:
        with open(j) as infile:
            data = json.load(infile)
        # get the values under the "queries" key, add them to the merged dict
        merged_queries.update(data["queries"])

    ############## write json #################

    of3_query_dict = {"queries": merged_queries}

    with open(outpath, "w") as outfile:
        json.dump(of3_query_dict, outfile)

    print(f"Written merged json to {outpath}")


# =============================================================================
def main() -> None:
    """Workflow:
    ---
    main
     ├── args
     └── merge_many_jsons
    """
    # get args
    args = parse_args()

    # get all input jsons
    filepaths = Path(args.indir).glob("*.json")
    jsons = sorted([f for f in filepaths if f.is_file()])

    # func
    merge_many_jsons(json_list=jsons, outpath=args.outpath, args=args)


# =============================================================================
if __name__ == "__main__":
    sys.exit(main())
