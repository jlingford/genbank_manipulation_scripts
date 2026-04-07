#!/usr/bin/env python3
"""
Convert protein fasta file to OpenFold3 formatted query .json file for inference

OF3 requires a .json input in this format:

    {
        "queries": {
            "query_1": { ... },
            "query_2": { ... }
        }
    }

Where each query is a dict[str, list[dict]]:

        "query_1": {
            "chains": [ { ... }, { ... } ],
        }

And where each chain is a dict[str, str|list[str]], like so:

         {
             "molecule_type": "protein",
             "chain_ids": "A",
             "description": "Optional metadata example",
             "sequence": "PVLSCGEWQCL",
             "use_msas": true,
             "use_main_msas": true,
             "use_paired_msas": true,
             "main_msa_file_paths": "/absolute/path/to/main_msas",
             "paired_msa_file_paths": "/absolute/path/to/paired_msas",
             "template_alignment_file_path": "/absolute/path/to/template_msa",
             "template_entry_chain_ids": [
                 "entry1_A",
                 "entry2_B",
                 "entry3_A",
             ],
         },

Only one query is written per json file, which is

Input:
    - protein fasta file (.faa | .fasta)
        - assumes that one fasta file = one OF3 query input. i.e., a fasta file of multiple sequences is a multimer of different chains
        - will not work with the ColabFold multimer format (i.e., multimer chains separated by ":")
Output:
    - OF3 query .json file
        - one query only; query named after the input fasta filename
        - .json file named after the input .faa file

\033[1m\033[31mTip:\033[0m
    Run in parallel with:

        `parallel faa_to_of3json.py -i {} ::: *.faa`

\033[1m\033[31mWARNING:\033[0m
    ...

"""
# TODO:
# - [ ] add gzip .faa file reading?
# - [ ] improve chain labels to be infinite?
# - [ ] parse for identical sequences and make chain_ids a list[str]
# - [ ] parse description for msa.a3m filepath info and add to msa field of dict

from Bio import SeqIO
from pathlib import Path
from typing import TextIO, NamedTuple
import argparse
import gzip
import sys
import json
from collections import defaultdict


# =============================================================================
# CLI args
# =============================================================================
class Args(NamedTuple):
    infile: Path
    outdir: Path


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        type=Path,
        metavar="FILE",
        required=True,
        help="Path to input [Required]",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        required=False,
        default=".",
        metavar="DIR",
        help="Output target directory [Optional][Default: cwd]",
    )

    args = Args(**vars(parser.parse_args()))

    return args


# =============================================================================
# Core func.
# =============================================================================
def faa_to_of3json(
    input_faa: Path,
    outdir: Path,
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
    outpath = Path(outdir) / f"{input_faa.stem}.json"
    outdir.mkdir(parents=True, exist_ok=True)

    ########### parse faa file ###########

    # TODO:
    uniqseqs = defaultdict(list)
    for rec in SeqIO.parse(input_faa, "fasta"):
        uniqseqs[str(rec.seq)].append(rec.id)

    CHAIN_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"

    chains_list = []

    # create dict for each chain and add to list
    for i, rec in enumerate(SeqIO.parse(input_faa, "fasta")):
        # set fresh dict each iteration
        chain_dict = {}

        # get vars
        chain_id = CHAIN_LABELS[i]
        desc = str(rec.description)
        seq = str(rec.seq)

        # update dict fields
        chain_dict.update({"molecule_type": "protein"})
        chain_dict.update({"chain_ids": chain_id})
        chain_dict.update({"description": desc})
        chain_dict.update({"sequence": seq})

        # add it to list
        chains_list.append(chain_dict)

    #############################################
    # write chain chains_list to query dict

    # make the query_id the .faa filename
    query_id = f"{input_faa.stem}"

    # build openfold3 query dict
    of3_query_dict = {
        "queries": {
            query_id: {"chains": chains_list},
        },
    }

    ############## write json #################

    with open(outpath, "w") as outfile:
        json.dump(of3_query_dict, outfile)


# =============================================================================
def main() -> None:
    """Workflow:
    ---
    main
     ├── args
     └── faa_to_of3json
    """
    # get args
    args = parse_args()

    # func
    faa_to_of3json(input_faa=args.infile, outdir=args.outdir, args=args)


# =============================================================================
if __name__ == "__main__":
    sys.exit(main())
