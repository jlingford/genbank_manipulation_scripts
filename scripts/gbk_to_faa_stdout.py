#!/usr/bin/env python3
"""
Print all protein fasta sequences from a Genbank file

Input:
    - path to a Genbank file

Output:
    - protein fasta printed to stdout

\033[1m\033[31mWARNING:\033[0m
    ...

"""
# TODO:
# - [ ] add ability to read gzipped genbank_files

from Bio import SeqIO
from concurrent.futures import ProcessPoolExecutor
from functools import partial
from pathlib import Path
import argparse
import gzip
import sys


# =================================================================
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--input",
        dest="input",
        type=Path,
        metavar="FILE",
        required=True,
        help="Path to input Genbank file (.gbk or .gbff) [Required]",
    )

    args = parser.parse_args()

    return args


# =================================================================
def gbk_to_faa(genbank_file: Path) -> None:
    """Reads Genbank file and prints protein fasta to stdout

    Args:
        genbank_file (Path): path to genbank_file

    Returns:
        None: Prints protein fasta to stdout. Header line includes description
    """

    for rec in SeqIO.parse(genbank_file, "genbank"):
        for feat in rec.features:
            if feat.type != "CDS":
                continue
            # extract faa info from qualifiers dict[str, list[str]]
            faa_id = feat.qualifiers.get("locus_tag", [0])[0]
            faa_desc = feat.qualifiers.get("product", [0])[0]
            faa_seq = feat.qualifiers.get("translation", [0])[0]
            # print to stdout
            print(f">{faa_id} {faa_desc}\n{faa_seq}")


# =================================================================
def main() -> None:
    args = parse_args()
    gbk_to_faa(genbank_file=args.input)


if __name__ == "__main__":
    sys.exit(main())
