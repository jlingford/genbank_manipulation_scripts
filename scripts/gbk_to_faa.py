#!/usr/bin/env python3
"""
Write out all protein fasta sequences from a Genbank file

Input:
    - path to a Genbank file

Output:
    - protein fasta written to file

\033[1m\033[31mTip:\033[0m
    Run in parallel with:

        `parallel gbk_to_faa.py -i {} ::: *.gbk`

"""
# TODO:
# - [x] add ability to read gzipped genbank_files

from Bio import SeqIO
from pathlib import Path
from typing import TextIO
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

    parser.add_argument(
        "-o",
        "--outdir",
        dest="outdir",
        type=Path,
        metavar="DIR",
        required=False,
        help="Path output directory for faa files [Default: same dir as input Genbank file]",
    )

    parser.add_argument(
        "--no_desc",
        dest="no_desc",
        action="store_true",
        help="Write out fasta file with just the accession ID in the header, but no description",
    )

    args = parser.parse_args()

    # create default outdir from input dir
    if args.outdir is None:
        input_parentdir = Path(args.input.parent)
        args.outdir = input_parentdir

    return args


# =============================================================
def open_gz(file: Path) -> TextIO:
    """Utility function: open file, even if it is gzipped"""
    if file.suffix == ".gz":
        return gzip.open(file, "rt")
    else:
        return open(file, "r")


# =================================================================
def gbk_to_faa(
    genbank_file: Path,
    outdir: Path,
    args: argparse.Namespace,
) -> None:
    """Reads Genbank file and writes to file

    Args:
        genbank_file (Path): path to genbank_file

    Returns:
        None: writes output .faa file that shares the same name as the input genbank_file
    """
    ############# make output file ###################
    outfaa = Path(outdir) / f"{args.input.stem}.faa"
    if outfaa.exists():
        outfaa.unlink()
    outfaa.parent.mkdir(parents=True, exist_ok=True)

    ############# read genbank and write faa ###################
    with (
        open_gz(file=genbank_file) as infile,
        open(file=outfaa, mode="w") as outfile,
    ):
        for rec in SeqIO.parse(infile, "genbank"):
            for feat in rec.features:
                if feat.type != "CDS":
                    continue

                # extract faa info from qualifiers dict[str, list[str]]
                faa_id = feat.qualifiers.get("locus_tag", [0])[0]
                faa_seq = feat.qualifiers.get("translation", [0])[0]
                # WARN: key for descriptions can vary between genbanks
                # faa_desc = feat.qualifiers.get("product", [0])[0]
                faa_desc = feat.qualifiers.get("annotation", [0])[0]

                if args.no_desc:
                    outfile.write(f">{faa_id}\n{faa_seq}\n")
                else:
                    outfile.write(f">{faa_id} {faa_desc}\n{faa_seq}\n")


# =================================================================
def main() -> None:
    args = parse_args()
    gbk_to_faa(genbank_file=args.input, outdir=args.outdir, args=args)


if __name__ == "__main__":
    sys.exit(main())
