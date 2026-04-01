#!/usr/bin/env python3
"""
Write out all DNA sequences from a Genbank file

Input:
    - path to a Genbank file

Output:
    - DNA fasta written to file

\033[1m\033[31mTip:\033[0m
    Run in parallel with:

        `parallel gbk_to_fna.py -i {} ::: *.gbk`

"""
# TODO:
# - [x] add ability to read gzipped genbank_files

from Bio import SeqIO
from pathlib import Path
from typing import TextIO, NamedTuple
import argparse
import gzip
import sys


# =================================================================
# CLI args
# =================================================================
class Args(NamedTuple):
    infile: Path
    outdir: Path
    add_desc: bool


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--infile",
        dest="infile",
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
        help="Path output directory for fna files [Default: same dir as input Genbank file]",
    )

    parser.add_argument(
        "--add_desc",
        dest="add_desc",
        action="store_true",
        help="Write out fasta file with the description in the header, after the accession ID",
    )

    args = Args(**vars(parser.parse_args()))

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
def gbk_to_fna(
    genbank_file: Path,
    outdir: Path,
    args: Args,
) -> None:
    """Reads Genbank file and writes to file

    Args:
        genbank_file (Path): path to genbank_file

    Returns:
        None: writes output .fna file that shares the same name as the input genbank_file
    """
    ############# make output file ###################
    outfna = Path(outdir) / f"{args.infile.stem}.fna"
    if outfna.exists():
        outfna.unlink()
    outfna.parent.mkdir(parents=True, exist_ok=True)

    ############# read genbank and write fna ###################

    with (
        open_gz(file=genbank_file) as infile,
        open(file=outfna, mode="w") as outfile,
    ):
        for rec in SeqIO.parse(infile, "genbank"):
            id = rec.id
            desc = rec.description
            dna = str(rec.seq)

            if args.add_desc:
                outfile.write(f">{id} {desc}\n{dna}\n")
            else:
                outfile.write(f">{id}\n{dna}\n")


# =================================================================
def main() -> None:
    args = parse_args()
    gbk_to_fna(genbank_file=args.infile, outdir=args.outdir, args=args)


if __name__ == "__main__":
    sys.exit(main())
