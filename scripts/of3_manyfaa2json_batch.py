#!/usr/bin/env python3
"""
Converts a directory of protein fasta files to a single OpenFold3 formatted query .json file for inference

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

\033[1m\033[31mWARNING:\033[0m
    - make sure each input .fasta file has a unique name, this will be the query ID in the .json outfile
    - make sure each directory of .fasta files has a unique name, as this dir name will become the name of the .json outfile
    ...

"""
# TODO:
# - [ ] add gzip .faa file reading?
# - [ ] improve chain labels to be infinite?
# - [x] parse for identical sequences and make chain_ids a list[str]
# - [x] parse description for msa.a3m filepath info and add to msa field of dict
# - [ ] improve outpath writing and defaults

from Bio import SeqIO
from pathlib import Path
from collections import defaultdict
from dataclasses import dataclass
from typing import TextIO, NamedTuple, Any
import argparse
import gzip
import sys
import json


# =============================================================================
# CHAIN_LABELS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ"
CHAIN_LABELS = [
    "A",
    "B",
    "C",
    "D",
    "E",
    "F",
    "G",
    "H",
    "I",
    "J",
    "K",
    "L",
    "M",
    "N",
    "O",
    "P",
    "Q",
    "R",
    "S",
    "T",
    "U",
    "V",
    "W",
    "X",
    "Y",
    "Z",
    "AA",
    "AB",
    "AC",
    "AD",
    "AE",
    "AF",
    "AG",
    "AH",
    "AI",
    "AJ",
    "AK",
    "AL",
    "AM",
    "AN",
    "AO",
    "AP",
    "AQ",
    "AR",
    "AS",
    "AT",
    "AU",
    "AV",
    "AW",
    "AX",
    "AY",
    "AZ",
]


# =============================================================================
# CLI args
# =============================================================================
# class Args(NamedTuple):
@dataclass()
class Args:
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
        metavar="DIR",
        help="Output target directory [Default: path/to/inputdir/name_of_inputdir.json]",
    )

    args = Args(**vars(parser.parse_args()))

    if not args.indir.is_dir():
        parser.error("Error: input must be a target directory of fasta files")

    # TODO: fix
    # create default outdir from input dir
    if args.outpath is None:
        input_parentdir = Path(args.indir.parent)
        args.outpath = input_parentdir

    return args


# =============================================================================
# Core func.
# =============================================================================
def faa_to_of3json(
    input_faa: Path,
    args: Args,
) -> dict[str, dict]:
    """Converts protein fasta info to json format for OpenFold3

    ---
    Args:
        input_faa (Path): path to .faa file
        outpath (Path): path to output dir
        args (Args): other args

    Returns:
        query_dict (dict): dict of the fasta query
    """
    ############# step 1 #################################

    # collect chain info, combine identical chain seqs on the same key
    uniqseqs = defaultdict(list)
    msapaths = {}
    metainfo = {}

    for i, rec in enumerate(SeqIO.parse(input_faa, "fasta")):
        chain_id = CHAIN_LABELS[i]
        desc = str(rec.description)
        seq = str(rec.seq)
        msa = str(rec.description.split("|")[-1])

        # update dicts
        uniqseqs[seq].append(chain_id)
        msapaths[seq] = msa
        metainfo[seq] = desc

    ############# step 2 #################################

    # unpack info from dicts above on shared sequence, add to new chain_dict
    # create dict for each chain and add to list
    chains_list = []
    for k in uniqseqs:
        # set fresh dict each iteration!
        chain_dict = {}

        # if only one chain id per chain, make it a str instead of list
        chain_id_format = uniqseqs[k][0] if len(uniqseqs[k]) == 1 else uniqseqs[k]

        # update dict fields
        chain_dict.update({"molecule_type": "protein"})
        chain_dict.update({"chain_ids": chain_id_format})
        chain_dict.update({"description": metainfo[k]})
        chain_dict.update({"main_msa_file_paths": msapaths[k]})
        chain_dict.update({"sequence": k})

        # add it to list
        chains_list.append(chain_dict)

    ############# step 3 #################################

    # make the query_id the .faa filename
    query_id = f"{input_faa.stem}"

    # construct dict for query id
    query_dict: dict[str, dict[str, list[Any]]] = {query_id: {"chains": chains_list}}

    return query_dict


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

    # get input fasta files in target dir
    faa_files = sorted([f for f in Path(args.indir).glob("*.fa*") if f.is_file()])

    # append each query dict on top of each other, where each key is the faa file name
    queries = {}
    for faa in faa_files:
        new_query = faa_to_of3json(input_faa=faa, args=args)
        queries.update(new_query)

    # combine all queries into the of3 dict
    all_queries = {"queries": queries}

    # set outpath
    of3_json_outpath = (
        Path(args.outpath)
        / f"{args.indir.stem.removesuffix('.json')}.json"  # remove any suffix incase user added one, and add it back
    )
    of3_json_outpath.parent.mkdir(parents=True, exist_ok=True)
    # write json
    with open(of3_json_outpath, "w") as of:
        json.dump(all_queries, of)

    print(f"OpenFold3 .json query file written to: {of3_json_outpath}")


# =============================================================================
if __name__ == "__main__":
    sys.exit(main())
