#!/usr/bin/env python3
"""
Plot jointplot with cutoff scores

Input:
    - table of genbank slices sizes and hyddb pident scores
        - .tsv format
        - with header

        col1: gene_id
        col2: hyddb_group name
        col3: pident
        col4: classification_confidence
        col5: size_bp

HARDCODED VARS:
    - cutoff scores x/y values
"""

import sys
import argparse
import numpy as np
import polars as pl
import seaborn as sns
from pathlib import Path
from typing import NamedTuple, TextIO
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
from matplotlib import rcParams


# =============================================================================
# GLOBAL VARS
# =============================================================================

plottype = "jointplot"

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
        default=f"{plottype}",
        metavar="PATH",
        help=f"Output target directory [Optional][Default: ./{plottype}/]",
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
        has_header=True,
    ).select(
        [
            # "gene_id",
            # "hyddb_group",
            "pident",
            # "classification_confidence",
            "size_bp",
        ]
    )
    print(df)

    ############### write data #####################
    name, _ = themeing(args=args)
    outpath = Path(outpath.parent) / f"{outpath.stem}-{name}-{plottype}"
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

    # xaxis from df
    # WARN: hardcoded
    x = df["size_bp"]
    y = df["pident"]

    # axis labels
    # WARN: hardcoded
    ylabel = "HydDB closest hit seq. id. [%]"
    xlabel = "contig length [bp]"

    # makes a dataframe of sample sizes. get the number out of the list with by calling index [0]
    size_n = len(df)

    # add labels with sample size
    label = f"{name}\n$\\mathit{{n}}$ = {size_n:,}"

    # size cutoff line to annotate on histplot
    xline_cutoff = 10000
    yline_cutoff = 70
    # if condition meets cutoff, then set alpha to 0.3, otherwise 0.9
    condition = (y > 70) | (x < 10000)
    alphas = np.where(condition, 0.30, 0.9)

    # make bin size and xmin/xmax
    # WARN: hardcoded
    BINS = 50
    xmin, xmax = 0, 25000
    bins = np.linspace(xmin, xmax, BINS + 1)

    ################## create plot ########################
    fig, ax = plt.subplots(
        nrows=1,
        ncols=1,
        figsize=(
            15,  # width
            15,  # height
        ),
    )

    jp = sns.jointplot(
        df,
        x=x,
        y=y,
        kind="scatter",
        # space=0,  # space between marginal histplot and scatterplot
        # alpha=0.8,
        # size=2,  # size of scatterpoints
        edgecolor=None,
        linewidth=0,
        color=color,
        # marginal_kws={"color": "#a6adc8"},  # marginal histplot settings
        # joint_kws=dict(
        #     gridsize=50,
        #     mincnt=1,
        # ),
        # cmap=cmap,
    )

    # custom x/y labels
    jp.set_axis_labels(
        ylabel=ylabel,
        xlabel=xlabel,
    )

    # call axes of jointplot "ax"
    ax = jp.ax_joint

    ################## scatter aesthetics ##################

    # change scatter colors based on cutoff conditions
    scatter = jp.ax_joint.collections[0]
    facecolors = scatter.get_facecolor()  # preserve the existing color

    # facecolors need to be interpreted as a single row, use np.tile
    facecolors = np.tile(facecolors, (len(x), 1))

    # update alphas of scatter points
    facecolors[:, 3] = alphas
    scatter.set_facecolor(facecolors)

    ################## plot aesthetics ####################

    # add sample size to top left corner
    ax.text(
        0.02,  # x
        0.98,  # y
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
    )

    # add cutoff lines
    ax.axvline(x=xline_cutoff, color="grey", linewidth=1, linestyle="--")
    ax.axhline(y=yline_cutoff, color="grey", linewidth=1, linestyle="--")

    ################## write plot ###################
    outpath = Path(outpath.parent) / f"{outpath.stem}-{name}-{plottype}"
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
