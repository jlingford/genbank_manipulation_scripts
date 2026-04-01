#!/usr/bin/env python3
"""
Plot histogram

Input:
    - table of genbank slices sizes (.tsv, no header)

        col1: gene_names
        col: size in bp (int)
"""

import sys
import argparse
import numpy as np
import polars as pl
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
from typing import NamedTuple, TextIO
from matplotlib import rcParams
import matplotlib.font_manager as fm


# =============================================================================
# GLOBAL VARS
# =============================================================================

# FONTS
arial_font = "../auxfiles/arial.ttf"
arial_font_bold = "../auxfiles/Arial Bold.ttf"
fm.fontManager.addfont(arial_font)
fm.fontManager.addfont(arial_font_bold)
rcParams["font.sans-serif"] = "Arial"
rcParams["font.family"] = "Arial"
rcParams["font.size"] = 10


# =============================================================================
# CLI args
# =============================================================================
class Args(NamedTuple):
    infile: Path
    outpath: Path
    dpi: int
    nife: bool
    fefe: bool


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
        help="Path to input file [REQUIRED]",
    )

    parser.add_argument(
        "-o",
        "--outpath",
        type=Path,
        required=False,
        default="plot_output",
        metavar="PATH",
        help="Output target directory [Optional][Default: ./plot_output]",
    )

    parser.add_argument(
        "-D",
        "--dpi",
        dest="dpi",
        type=int,
        required=False,
        default=300,
        metavar="N",
        help="Resolution in dpi for output .png [Default: 300]",
    )

    parser.add_argument(
        "-N",
        dest="nife",
        action="store_true",
        required=False,
        help="Set NiFe name and preferred colorscheme (convenience setting to avoid manually editing this script)",
    )

    parser.add_argument(
        "-F",
        dest="fefe",
        action="store_true",
        required=False,
        help="Set FeFe name and preferred colorscheme (convenience setting to avoid manually editing this script)",
    )

    args = Args(**vars(parser.parse_args()))

    return args


# =============================================================================
def themeing(args: Args) -> tuple[str, str]:
    """Util function: get plot name and colorscheme"""
    # default
    name = f"{args.infile.stem}"
    color = "#b4befe"
    # setting from CLI args
    if args.nife:
        name = "NiFe"
        color = "#94e2d5"
    if args.fefe:
        name = "FeFe"
        color = "#fab387"

    return name, color


# =============================================================================
def process_data(infile: Path, outpath: Path, args: Args) -> pl.Dataframe:
    """Process raw data and return polars dataframe"""
    ############### create outpath ################
    outpath.parent.mkdir(parents=True, exist_ok=True)

    ############### read data #####################
    df = pl.read_csv(
        infile,
        separator="\t",
        has_header=False,
        new_columns=[
            "gene_name",
            "sequence_length",
        ],
    ).select(
        [
            "sequence_length",
        ]
    )
    print(df)

    ############### write data #####################
    name, _ = themeing(args=args)
    outpath = Path(outpath.parent) / f"{outpath.stem}-{name}-histplot"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    df.write_csv(
        f"{outpath}_dataframe.tsv",
        separator="\t",
        include_header=True,
    )

    return df


# =============================================================================
def plot(
    df: pl.DataFrame,
    outpath: Path,
    args: Args,
) -> None:
    """Make histplot"""
    ############ Prepare figure elements ############

    # get name and color
    name, color = themeing(args=args)

    # axis labels
    # WARN: hardcoded
    ylabel = "count"
    xlabel = "sequence length [bp]"

    # makes a dataframe of sample sizes. get the number out of the list with by calling index [0]
    size_n = len(df)

    # add labels with sample size
    label = f"{name}\n$\\mathit{{n}}$ = {size_n:,}"

    # xaxis from df
    # WARN: hardcoded
    x = "sequence_length"
    # y = "sequence_length"

    # make bin size and xmin/xmax
    # WARN: hardcoded
    BINS = 50
    xmin, xmax = 0, 25000
    bins = np.linspace(xmin, xmax, BINS + 1)

    # size cutoff line to annotate on histplot
    vline_cutoff = 10000

    ################## create plot ########################
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(
            10,  # width
            5,  # height
        ),
    )

    sns.histplot(
        df,
        x=x,
        # y=y,
        ax=ax,  # use the seaborn ax parameter to make the subplots work
        bins=bins,
        color=color,
        zorder=3,
        # alpha=1,
        # NOTE: remove lines of histogram bars and add them in later to avoid the artifact line the bins array introduces
        edgecolor=None,
        linewidth=0,
    )

    ################## plot aesthetics ####################

    # change line color and width on each histogram bar
    for bar in ax.patches:
        if bar.get_height() > 0:
            bar.set_edgecolor("black")
            bar.set_linewidth(0.5)

    # add sample size to top left corner
    ax.text(
        0.02,  # x
        1.10,  # y
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
    )

    # add cutoff line
    ax.axvline(x=vline_cutoff, color="grey", linewidth=1, linestyle="--")

    # spines and labels
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)
    ax.spines["left"].set_visible(True)
    ax.set_ylim(0, top=None)
    ax.set_xlim(xmin, xmax)
    # ax.xaxis.set_tick_params(length=0)  # remove ticks but keep labels
    ax.set_ylabel(ylabel)
    ax.set_xlabel(xlabel)
    # ax.grid(axis="x", color="grey", linestyle="-", linewidth=0.5, zorder=0)
    # ax.spines["bottom"].set_visible(False)
    # ax.tick_params(bottom=False, labelbottom=False)
    ax.spines["bottom"].set_visible(True)

    ################## write plot ###################
    outpath = Path(outpath.parent) / f"{outpath.stem}-{name}-histplot"
    outpath.parent.mkdir(parents=True, exist_ok=True)

    plt.savefig(
        f"{outpath}.png",
        format="png",
        dpi=args.dpi,
        bbox_inches="tight",
    )
    plt.savefig(
        f"{outpath}.svg",
        format="svg",
        bbox_inches="tight",
    )
    plt.savefig(
        f"{outpath}.pdf",
        format="pdf",
        bbox_inches="tight",
    )

    ################## show plot ###################
    plt.tight_layout()
    plt.show()
    plt.close("all")


# =============================================================================
def main() -> None:
    # get args
    args = parse_args()

    # create dataframe
    df = process_data(infile=args.infile, outpath=args.outpath, args=args)

    # generate plot
    plot(df=df, outpath=args.outpath, args=args)


# =============================================================================
if __name__ == "__main__":
    sys.exit(main())
