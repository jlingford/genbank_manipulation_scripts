#!/usr/bin/env python3
"""
BATCH writes out all DNA fasta sequences from a Genbank file

Input:
    - path to input directory containing Genbank files

Output:
    - DNA fasta (.fna)
"""
# TODO:
# - [x] add ability to read gzipped genbank_files
# - [ ] add regex to filter input files to .gbk or .gbff

from concurrent.futures import ProcessPoolExecutor
from functools import partial
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
    indir: Path
    outdir: Path
    cpu: int
    add_desc: bool


def parse_args() -> Args:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "-i",
        "--indir",
        type=Path,
        metavar="DIR",
        required=True,
        help="Path to target directory containing all input Genbank files for conversion (.gbk or .gbff) [Required]",
    )

    parser.add_argument(
        "-o",
        "--outdir",
        type=Path,
        metavar="DIR",
        required=False,
        help="Path output directory for fna files [Default: same dir as input dir]",
    )

    parser.add_argument(
        "-c",
        "--cpu",
        type=int,
        default=None,
        metavar="N",
        required=False,
        help="No. of CPUs to use for parallelism [Default: max available]",
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
        args.outdir = Path(args.input_dir)

    return args


# =============================================================
# Util
# =============================================================
def open_gz(file: Path) -> TextIO:
    """Utility function: open file, even if it is gzipped"""
    if file.suffix == ".gz":
        return gzip.open(file, "rt")
    else:
        return open(file, "r")


# =============================================================
# Core functions
# =================================================================
def gbk_to_fna(
    genbank_file: Path,
    outdir: Path,
    args: Args,
) -> None:
    """Reads Genbank file and writes out DNA to fasta file

    ---
    Args:
        genbank_file (Path): path to genbank_file

    Returns:
        None: writes output .fna file that shares the same name as the input genbank_file
    """
    ############# make output file ###################
    outfna = Path(outdir) / f"{genbank_file.stem}.fna"
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
    # parse args
    args = parse_args()

    # get all input genbank files
    gbk_files = sorted([f for f in Path(args.indir).glob("*.gb*") if f.is_file()])

    ############### PARALLEL PROCESSING ###################

    # gbk_to_fna(genbank_file=args.input, outdir=args.outdir, args=args)

    # make partial func
    partial_gbk_to_fna = partial(
        gbk_to_fna,
        outdir=args.outdir,
        args=args,
    )
    # ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=args.cpu) as exe:
        list(exe.map(partial_gbk_to_fna, gbk_files))


# =================================================================
if __name__ == "__main__":
    sys.exit(main())
